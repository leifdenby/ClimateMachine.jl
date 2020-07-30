# # Heat equation tutorial - Compare time steppers

# In this tutorial, we'll be solving the [heat
# equation](https://en.wikipedia.org/wiki/Heat_equation):

# ``
# \frac{∂ ρcT}{∂ t} + ∇ ⋅ (-α ∇ρcT) = 0
# ``

# where
#  - `t` is time
#  - `α` is the thermal diffusivity
#  - `T` is the temperature
#  - `ρ` is the density
#  - `c` is the heat capacity
#  - `ρcT` is the thermal energy

# To put this in the form of ClimateMachine's [`BalanceLaw`](@ref
# ClimateMachine.BalanceLaws.BalanceLaw), we'll re-write the equation as:

# ``
# \frac{∂ ρcT}{∂ t} + ∇ ⋅ (F(α, ρcT, t)) = 0
# ``

# where
#  - ``F(α, ρcT, t) = -α ∇ρcT`` is the second-order flux

# with boundary conditions
#  - Fixed temperature ``T_{surface}`` at ``z_{min}`` (non-zero Dirichlet)
#  - No thermal flux at ``z_{min}`` (zero Neumann)

# Solving these equations is broken down into the following steps:
# 1) Preliminary configuration
# 2) PDEs
# 3) Space discretization
# 4) Time discretization
# 5) Solver hooks / callbacks
# 6) Solve
# 7) Post-processing

# # Preliminary configuration

# ## [Loading code](@id Loading-code-heat)

# First, we'll load our pre-requisites:
#  - load external packages:
using MPI
using OrderedCollections
using Plots
using StaticArrays
using OrdinaryDiffEq
using DiffEqBase

#  - load CLIMAParameters and set up to use it:

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

#  - load necessary ClimateMachine modules:
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.BalanceLaws:
    BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux

using ClimateMachine.Mesh.Geometry: LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils

#  - import necessary ClimateMachine modules: (`import`ing enables us to
#  provide implementations of these structs/methods)
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

# ## Initialization

# Define the float type (`Float64` or `Float32`)
FT = Float64;
# Initialize ClimateMachine for CPU.
ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));

# Load some helper functions for plotting
include(joinpath(clima_dir, "docs", "plothelpers.jl"));

# # Define the set of Partial Differential Equations (PDEs)

# ## Define the model

# Model parameters can be stored in the particular [`BalanceLaw`](@ref
# ClimateMachine.BalanceLaws.BalanceLaw), in this case, a `HeatModel`:

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

# This model dictates the flow control, using [Dynamic Multiple
# Dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch), for which
# kernels are executed.

# ## Define the variables

# All of the methods defined in this section were `import`ed in
# [Loading code](@ref Loading-code-heat) to let us provide
# implementations for our `HeatModel` as they will be used by
# the solver.

# Specify auxiliary variables for `HeatModel`
vars_state(::HeatModel, ::Auxiliary, FT) = @vars(z::FT, T::FT);

# Specify prognostic variables, the variables solved for in the PDEs, for
# `HeatModel`
vars_state(::HeatModel, ::Prognostic, FT) = @vars(ρcT::FT);

# Specify state variables whose gradients are needed for `HeatModel`
vars_state(::HeatModel, ::Gradient, FT) = @vars(ρcT::FT);

# Specify gradient variables for `HeatModel`
vars_state(::HeatModel, ::GradientFlux, FT) = @vars(α∇ρcT::SVector{3, FT});

# ## Define the compute kernels

# Specify the initial values in `aux::Vars`, which are available in
# `init_state_prognostic!`. Note that
# - this method is only called at `t=0`
# - `aux.z` and `aux.T` are available here because we've specified `z` and `T`
# in `vars_state` given `Auxiliary`
function init_state_auxiliary!(m::HeatModel, aux::Vars, geom::LocalGeometry)
    aux.z = geom.coord[3]
    aux.T = m.initialT
end;

# Specify the initial values in `state::Vars`. Note that
# - this method is only called at `t=0`
# - `state.ρcT` is available here because we've specified `ρcT` in
# `vars_state` given `Prognostic`
function init_state_prognostic!(
    m::HeatModel,
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
)
    state.ρcT = m.ρc * aux.T
end;

# The remaining methods, defined in this section, are called at every
# time-step in the solver by the [`BalanceLaw`](@ref
# ClimateMachine.BalanceLaws.BalanceLaw) framework.

# Overload `update_auxiliary_state!` to call `heat_eq_nodal_update_aux!`, or
# any other auxiliary methods
function update_auxiliary_state!(
    dg::DGModel,
    m::HeatModel,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
    nodal_update_auxiliary_state!(heat_eq_nodal_update_aux!, dg, m, Q, t, elems)
end;

# Compute/update all auxiliary variables at each node. Note that
# - `aux.T` is available here because we've specified `T` in
# `vars_state` given `Auxiliary`
function heat_eq_nodal_update_aux!(
    m::HeatModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    aux.T = state.ρcT / m.ρc
end;

# Since we have second-order fluxes, we must tell `ClimateMachine` to compute
# the gradient of `ρcT`. Here, we specify how `ρcT` is computed. Note that
#  - `transform.ρcT` is available here because we've specified `ρcT` in
#  `vars_state` given `Gradient`
function compute_gradient_argument!(
    m::HeatModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    transform.ρcT = state.ρcT
end;

# Specify where in `diffusive::Vars` to store the computed gradient from
# `compute_gradient_argument!`. Note that:
#  - `diffusive.α∇ρcT` is available here because we've specified `α∇ρcT` in
#  `vars_state` given `Gradient`
#  - `∇transform.ρcT` is available here because we've specified `ρcT`  in
#  `vars_state` given `Gradient`
function compute_gradient_flux!(
    m::HeatModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    diffusive.α∇ρcT = -m.α * ∇transform.ρcT
end;

# We have no sources, nor non-diffusive fluxes.
function source!(m::HeatModel, _...) end;
function flux_first_order!(m::HeatModel, _...) end;

# Compute diffusive flux (``F(α, ρcT, t) = -α ∇ρcT`` in the original PDE).
# Note that:
# - `diffusive.α∇ρcT` is available here because we've specified `α∇ρcT` in
# `vars_state` given `GradientFlux`
function flux_second_order!(
    m::HeatModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
    flux.ρcT += diffusive.α∇ρcT
end;

# ### Boundary conditions

# Second-order terms in our equations, ``∇⋅(F)`` where ``F = -α∇ρcT``, are
# internally reformulated to first-order unknowns.
# Boundary conditions must be specified for all unknowns, both first-order and
# second-order unknowns which have been reformulated.

# The boundary conditions for `ρcT` (first order unknown)
function boundary_state!(
    nf,
    m::HeatModel,
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
    ## Apply Dirichlet BCs
    if bctype == 1 # At bottom
        state⁺.ρcT = m.ρc * m.T_bottom
    elseif bctype == 2 # At top
        nothing
    end
end;

# The boundary conditions for `ρcT` are specified here for second-order
# unknowns
function boundary_state!(
    nf,
    m::HeatModel,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
    ## Apply Neumann BCs
    if bctype == 1 # At bottom
        nothing
    elseif bctype == 2 # At top
        diff⁺.α∇ρcT = n⁻ * m.flux_top
    end
end;

# # Spatial discretization

# Prescribe polynomial order of basis functions in finite elements
N_poly = 5;

# Specify the number of vertical elements
nelem_vert = 10;

# Specify the domain height
zmax = FT(1);

# Specify simulation time (SI units)
t0 = FT(0)
timeend = FT(2)

# Here, we define a convenience function to
# collect the prognostic and auxiliary states
# from the solver configuration.
function dict_of_states(solver_config, z_key)
    FT = eltype(solver_config.Q)
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        solver_config.dg.grid,
        solver_config.Q,
        vars_state(solver_config.dg.balance_law, Prognostic(), FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        solver_config.dg.grid,
        solver_config.dg.state_auxiliary,
        vars_state(solver_config.dg.balance_law, Auxiliary(), FT);
        exclude = [z_key],
    )
    return OrderedDict(state_vars..., aux_vars...)
end

z_scale = 100; # convert from meters to cm
z_key = "z";
z_label = "z [cm]";

output_dir = @__DIR__;
mkpath(output_dir);

# Define the number of outputs from `t0` to `timeend`
const n_outputs = 5;

# This equates to exports every ceil(Int, timeend/n_outputs) time-step:
const every_x_simulation_time = ceil(Int, timeend / n_outputs);

z_coord = FT[]
perf = Dict()
α_all = [2^FT(k) for k in -8:-5]
solver_types = [ImplicitSolverType(OrdinaryDiffEq.ROCK4()), ExplicitSolverType()]

all_data = Dict()  # store initial condition at ``t=0``
time_data = Dict()  # store time data

for α in α_all

    # Create an instance of the `HeatModel`:
    m = HeatModel{FT, typeof(param_set)}(; param_set = param_set, α=α);

    # Establish a `ClimateMachine` single stack configuration
    driver_config = ClimateMachine.SingleStackConfiguration(
        "HeatEquation",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
    );

    Δ = min_node_distance(driver_config.grid)
    given_Fourier = FT(0.08);
    Fourier_bound = given_Fourier * Δ^2 / m.α;
    dt = Fourier_bound
    z_coord = get_z(driver_config.grid, z_scale);

    for solver_type in solver_types
        case = (typeof(solver_type), α)
        @info case

        solver_config = ClimateMachine.SolverConfiguration(
                t0,
                timeend,
                driver_config,
                ode_dt = dt;
                ode_solver_type = solver_type)

        all_data[case...] = Dict[dict_of_states(solver_config, z_key)]  # store initial condition at ``t=0``
        time_data[case...] = FT[0]                                      # store time data

        callback = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
            push!(all_data[case...], dict_of_states(solver_config, z_key))
            push!(time_data[case...], gettime(solver_config.solver))
            nothing
        end;

        perf[case...] =
            @elapsed ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));

        # Append result at the end of the last time step:
        push!(all_data[case...], dict_of_states(solver_config, z_key))
        push!(time_data[case...], gettime(solver_config.solver))
    end # time-steppers
end # thermal diffusivity (α)


scaling = [perf[typeof(solver_types[1]), α] for α in α_all]
plot(α_all, scaling, label="$(typeof(solver_types[1]))")
scaling = [perf[typeof(solver_types[2]), α] for α in α_all]
plot!(α_all, scaling, label="$(typeof(solver_types[2]))")
plot!(xlabel="α (stiffness)", ylabel="runtime", title="Stepping comparison")
savefig("performance_comparison.png")
# ![](performance_comparison.png)

plot()
for α in α_all
    for solver_type in solver_types
        case = (typeof(solver_type), α)
        plot!(all_data[case...][end]["T"], z_coord, label="case=$case, T")
    end
end
savefig("steady_state_solutions.png")
# ![](steady_state_solutions.png)


# ```julia
# include(joinpath("tutorials", "Land", "Heat", "compare_steppers.jl"))
# ```
