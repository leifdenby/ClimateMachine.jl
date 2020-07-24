using ClimateMachine
using ClimateMachine.SystemSolvers
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
import ClimateMachine.ODESolvers:gettime, getdt
using ClimateMachine.SingleStackUtils
using ClimateMachine.BalanceLaws:
    BalanceLaw,
    Auxiliary,
    Gradient,
    GradientFlux,
    Prognostic
import ClimateMachine.BalanceLaws:
    vars_state,
    source!,
    flux_second_order!,
    flux_first_order!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    update_auxiliary_state!,
    nodal_update_auxiliary_state!,
    init_state_auxiliary!,
    init_state_prognostic!,
    boundary_state!

using MPI
using Plots
using OrderedCollections
using StaticArrays
using DiffEqBase
using OrdinaryDiffEq: Rosenbrock23
using OrdinaryDiffEq
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

FT = Float64;

ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));

include(joinpath(clima_dir, "docs", "plothelpers.jl"));

Base.@kwdef struct HeatModel{FT, APS} <: BalanceLaw
    "Parameters"
    param_set::APS
    "Heat capacity"
    ρc::FT = 1
    "Thermal diffusivity"
    α::FT = 0.01
    "Initial conditions for temperature"
    initialT::FT = 295.15
    "Bottom boundary value for temperature (Dirichlet boundary conditions)"
    T_bottom::FT = 300.0
    "Top flux (α∇ρcT) at top boundary (Neumann boundary conditions)"
    flux_top::FT = 0.0
end

m = HeatModel{FT, typeof(param_set)}(;param_set=param_set);

vars_state(::HeatModel, ::Auxiliary, FT) = @vars(z::FT, T::FT);
vars_state(::HeatModel, ::Prognostic, FT) = @vars(ρcT::FT);
vars_state(::HeatModel, ::Gradient, FT) = @vars(ρcT::FT);
vars_state(::HeatModel, ::GradientFlux, FT) = @vars(α∇ρcT::SVector{3, FT});

function init_state_auxiliary!(m::HeatModel, aux::Vars, geom::LocalGeometry)
    aux.z = geom.coord[3]
    aux.T = m.initialT
end;
function init_state_prognostic!(m::HeatModel,state::Vars,aux::Vars,coords,t::Real)
    state.ρcT = m.ρc * aux.T
end;
function update_auxiliary_state!(dg::DGModel,m::HeatModel,Q,t::Real,elems::UnitRange)
    nodal_update_auxiliary_state!(heat_eq_nodal_update_aux!, dg, m, Q, t, elems)
end;
function heat_eq_nodal_update_aux!(m::HeatModel,state::Vars,aux::Vars,t::Real)
    aux.T = state.ρcT / m.ρc
end;
function compute_gradient_argument!(m::HeatModel,transform::Vars,state::Vars,aux::Vars,t::Real)
    transform.ρcT = state.ρcT
end;
function compute_gradient_flux!(m::HeatModel,diffusive::Vars,∇transform::Grad,state::Vars,aux::Vars,t::Real)
    diffusive.α∇ρcT = -m.α * ∇transform.ρcT
end;
function source!(m::HeatModel, _...) end;
function flux_first_order!(m::HeatModel,flux::Grad,state::Vars,aux::Vars,t::Real,direction) end;
function flux_second_order!(m::HeatModel,flux::Grad,state::Vars,diffusive::Vars,hyperdiffusive::Vars,aux::Vars,t::Real)
    flux.ρcT += diffusive.α∇ρcT
end;
function boundary_state!(nf,m::HeatModel,state⁺::Vars,aux⁺::Vars,n⁻,state⁻::Vars,aux⁻::Vars,bctype,t,_...)
    bctype == 1 && (state⁺.ρcT = m.ρc * m.T_bottom)
end;
function boundary_state!(nf,m::HeatModel,state⁺::Vars,diff⁺::Vars,aux⁺::Vars,n⁻,state⁻::Vars,diff⁻::Vars,aux⁻::Vars,bctype,t,_...)
    bctype == 2 && (diff⁺.α∇ρcT = n⁻ * m.flux_top)
end;

n_poly = 5;
nelem_vert = 10;
zmax = FT(1);

driver_config = ClimateMachine.SingleStackConfiguration("HeatEquation",
    n_poly, nelem_vert, zmax, param_set, m, numerical_flux_first_order = CentralNumericalFluxFirstOrder());

t0 = FT(0)
timeend = FT(40)
Δ = min_node_distance(driver_config.grid)
given_Fourier = FT(0.08);
Fourier_bound = given_Fourier * Δ^2 / m.α;
dt = Fourier_bound

sc = ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);

prob = ODEProblem(sc.dg, sc.Q, (0.0, sc.timeend),nothing)

mutable struct DiffEqJLSolver{I} <: ODESolvers.AbstractDiffEqJLSolver
    integ::I
end
integrator = DiffEqBase.init(prob,
    Kvaerno3(autodiff=false,linsolve=LinSolveGMRES()),
    save_everystep = false,
    save_start = false,
    save_end = false,)
solver = DiffEqJLSolver(integrator)
gettime(solver::DiffEqJLSolver) = solver.integ.t
getdt(solver::DiffEqJLSolver) = solver.integ.dt

grid = sc.dg.grid;
Q = sc.Q;
aux = sc.dg.state_auxiliary;

output_dir = @__DIR__;

mkpath(output_dir);

z_scale = 100 # convert from meters to cm
z_key = "z"
z_label = "z [cm]"
z = get_z(grid, z_scale)
st_prog = vars_state(m, Prognostic(), FT)
st_aux = vars_state(m, Auxiliary(), FT)
state_vars = get_vars_from_nodal_stack(grid,Q,st_prog)
aux_vars = get_vars_from_nodal_stack(grid,aux,st_aux)
all_vars = OrderedDict(state_vars..., aux_vars...);
f = joinpath(output_dir, "initial_condition.png")
export_plot_snapshot(z, all_vars, ("ρcT",), f, z_label);

n_outputs = 5;
every_x_simulation_time = ceil(Int, timeend / n_outputs);
all_data = Dict[Dict([k => Dict() for k in 0:n_outputs]...),]
all_data[1] = all_vars # store initial condition at ``t=0``

callback = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        grid,
        Q,
        st_prog,
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        grid,
        aux,
        st_aux;
        exclude = [z_key],
    )
    all_vars = OrderedDict(state_vars..., aux_vars...)
    push!(all_data, all_vars)

    nothing
end;

@time ODESolvers.solve!(sc.Q, solver; timeend = sc.timeend, callbacks = (callback,))

@show keys(all_data[1])

f = joinpath(output_dir, "solution_vs_time.png");
export_plot(z, all_data, ("ρcT",), f, z_label);
