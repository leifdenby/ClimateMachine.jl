# # Heat equation tutorial

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
using Statistics
using Dierckx

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

# Create an instance of the `HeatModel`:
m = HeatModel{FT, typeof(param_set)}(; param_set = param_set);

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
        nothing#state⁺.ρcT = m.ρc * m.T_bottom
    elseif bctype == 2 # At top
        state⁺.ρcT = m.ρc * m.T_bottom#nothing
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
        diff⁺.α∇ρcT = n⁻ * m.flux_top #nothing
    elseif bctype == 2 # At top
        nothing #diff⁺.α∇ρcT = n⁻ * m.flux_top
    end
end;

# # Spatial discretization

# Prescribe polynomial order of basis functions in finite elements
N_poly = 5;

# Specify the number of vertical elements
nelem_vert = 10;

# # Specify the domain height
# zmax = FT(1);
# Specify the domain boundaries
    zmax = FT(0)
    zmin = FT(-1)

# Establish a `ClimateMachine` single stack configuration
driver_config = ClimateMachine.SingleStackConfiguration(
    "HeatEquation",
    N_poly,
    nelem_vert,
    zmax,
    param_set,
    m,
    zmin = zmin,
    numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
);

# # Time discretization

# Specify simulation time (SI units)
t0 = FT(0)
timeend = FT(1)

# We'll define the time-step based on the [Fourier
# number](https://en.wikipedia.org/wiki/Fourier_number)
Δ = min_node_distance(driver_config.grid)

given_Fourier = FT(0.08);
Fourier_bound = given_Fourier * Δ^2 / m.α;
dt = Fourier_bound

# # Configure a solver.

# This initializes the state vector and allocates memory for the solution in
# space (`dg` has the model `m`, which describes the PDEs as well as the
# function used for initialization). This additionally initializes the ODE
# solver, by default an explicit Low-Storage
# [Runge-Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
# method.

solver_config =
    ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
grid = solver_config.dg.grid;
Q = solver_config.Q;
aux = solver_config.dg.state_auxiliary;

# ## Inspect the initial conditions

# Let's export a plot of the initial state
output_dir = @__DIR__;

mkpath(output_dir);

z_scale = 100; # convert from meters to cm
z_key = "z";
z_label = "z [cm]";
z = get_z(grid, z_scale);

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

# Create an array to store the solution:
all_data = Dict[dict_of_states(solver_config, z_key)]  # store initial condition at ``t=0``
time_data = FT[0]                                      # store time data

export_plot(
    z,
    all_data,
    ("ρcT",),
    joinpath(output_dir, "initial_condition.png");
    xlabel = "ρcT",
    ylabel = z_label,
    time_data = time_data,
);
# ![](initial_condition.png)

# It matches what we have in `init_state_prognostic!(m::HeatModel, ...)`, so
# let's continue.

# # Solver hooks / callbacks

# Define the number of outputs from `t0` to `timeend`
const n_outputs = 5;

# This equates to exports every ceil(Int, timeend/n_outputs) time-step:
const every_x_simulation_time = ceil(Int, timeend / n_outputs);

# The `ClimateMachine`'s time-steppers provide hooks, or callbacks, which
# allow users to inject code to be executed at specified intervals. In this
# callback, a dictionary of prognostic and auxiliary states are appended to
# `all_data` for time the callback is executed. In addition, time is collected
# and appended to `time_data`.
callback = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
    push!(all_data, dict_of_states(solver_config, z_key))
    push!(time_data, gettime(solver_config.solver))
    nothing
end;

# # Solve

# This is the main `ClimateMachine` solver invocation. While users do not have
# access to the time-stepping loop, code may be injected via `user_callbacks`,
# which is a `Tuple` of callbacks in [`GenericCallbacks`](@ref ClimateMachine.GenericCallbacks).
 @time ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));

# Append result at the end of the last time step:
push!(all_data, dict_of_states(solver_config, z_key));
push!(time_data, gettime(solver_config.solver));

# # Post-processing

# Our solution is stored in the array of dictionaries `all_data` whose keys are
# the output interval. The next level keys are the variable names, and the
# values are the values along the grid:

# To get `T` at ``t=0``, we can use `T_at_t_0 = all_data[1]["T"][:]`
@show keys(all_data[1])

# Let's plot the solution:

# #40 sec bonan run
# #include("./helperfunc.jl")
# bonan_z = reverse([
#        0.000,
#       -0.847,
#       -2.542,
#       -4.237,
#       -5.932,
#       -7.627,
#       -9.322,
#      -11.017,
#      -12.712,
#      -14.407,
#      -16.102,
#      -17.797,
#      -19.492,
#      -21.186,
#      -22.881,
#      -24.576,
#      -26.271,
#      -27.966,
#      -29.661,
#      -31.356,
#      -33.051,
#      -34.746,
#      -36.441,
#      -38.136,
#      -39.831,
#      -41.525,
#      -43.220,
#      -44.915,
#      -46.610,
#      -48.305,
#      -50.000,
#      -51.695,
#      -53.390,
#      -55.085,
#      -56.780,
#      -58.475,
#      -60.169,
#      -61.864,
#      -63.559,
#      -65.254,
#      -66.949,
#      -68.644,
#      -70.339,
#      -72.034,
#      -73.729,
#      -75.424,
#      -77.119,
#      -78.814,
#      -80.508,
#      -82.203,
#      -83.898,
#      -85.593,
#      -87.288,
#      -88.983,
#      -90.678,
#      -92.373,
#      -94.068,
#      -95.763,
#      -97.458,
#      -99.153
#      ])

#     bonan_temperature = reverse([
#      300.000,
#      299.969,
#      299.908,
#      299.847,
#      299.786,
#      299.725,
#      299.664,
#      299.604,
#      299.543,
#      299.483,
#      299.424,
#      299.365,
#      299.306,
#      299.248,
#      299.190,
#      299.133,
#      299.077,
#      299.021,
#      298.966,
#      298.911,
#      298.858,
#      298.805,
#      298.753,
#      298.702,
#      298.652,
#      298.603,
#      298.554,
#      298.507,
#      298.461,
#      298.416,
#      298.372,
#      298.329,
#      298.288,
#      298.248,
#      298.208,
#      298.171,
#      298.134,
#      298.099,
#      298.065,
#      298.033,
#      298.002,
#      297.972,
#      297.944,
#      297.917,
#      297.892,
#      297.868,
#      297.845,
#      297.825,
#      297.805,
#      297.788,
#      297.772,
#      297.757,
#      297.744,
#      297.733,
#      297.723,
#      297.715,
#      297.708,
#      297.704,
#      297.700,
#      297.699
# ])

#     bonan_z = bonan_z ./ 100.0
#     # Create an interpolation from the Bonan data
#     bonan_temperature_continuous = Spline1D(bonan_z, bonan_temperature)
#     bonan_at_clima_z = [bonan_temperature_continuous(i) for i in z]

   # plot([all_data[end]["T"]], z, label = ["Clima Heat Tutorial"])
   # plot([all_vars["soil.heat.T"] all_data[end]["T"]], z, label = ["Clima Land/Heat Model" "Clima Heat Tutorial"])
    # xlabel!("Temperature [K]")
    # ylabel!("Depth [cm]")

    # all_vars["soil.heat.T"] 

export_plot(
    z,
    all_data ,
    ("ρcT",),
    joinpath(output_dir, "solution_vs_time.png");
    xlabel = "ρcT",
    ylabel = z_label,
    time_data = time_data,
);
#![](solution_vs_time.png)

# The results look as we would expect: a fixed temperature at the bottom is
# resulting in heat flux that propagates up the domain. To run this file, and
# inspect the solution in `all_data`, include this tutorial in the Julia REPL
# with:

# ```julia
# include(joinpath("tutorials", "Land", "Heat", "heat_equation.jl"))
# ```
