using MPI
using OrderedCollections
using StaticArrays
<<<<<<< HEAD
=======
using Test
>>>>>>> master

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land
<<<<<<< HEAD
using ClimateMachine.Land.SoilWaterParameterizations
=======
>>>>>>> master
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils

FT = Float64;
<<<<<<< HEAD

# When using ClimateMachine, we need to define a function that sets the initial
# state of our model run.

function init_soil_water!(land, state, aux, coordinates, time)
    FT = eltype(state)
    state.soil.water.ϑ = FT(land.soil.water.initialϑ(aux))
    state.soil.water.θ_ice = FT(land.soil.water.initialθ_ice(aux))
end;




=======
>>>>>>> master
ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));

struct HeatModel end

<<<<<<< HEAD
SoilParams = SoilParamSet(porosity = 0.75, Ksat = 0.0, S_s = 1e-3)
surface_state = (aux, t) -> FT(0.1)
bottom_state = (aux, t) -> FT(0.4)
ϑ_0 = (aux) -> abs(aux.z) * 0.3 + 0.1
soil_water_model = SoilWaterModel(
    FT;
    params = SoilParams,
    initialϑ = ϑ_0,
    dirichlet_bc = Dirichlet(
        surface_state = surface_state,
        bottom_state = bottom_state,
    ),
    neumann_bc = Neumann(surface_flux = nothing, bottom_flux = nothing),
)

m_soil = SoilModel(soil_water_model, HeatModel())
sources = ()
m = LandModel(
    param_set,
    m_soil;
    source = sources,
    init_state_conservative = init_soil_water!,
)


N_poly = 5;
nelem_vert = 10;


# Specify the domain boundaries
zmax = FT(0);
zmin = FT(-1)

driver_config = ClimateMachine.SingleStackConfiguration(
    "LandModel",
    N_poly,
    nelem_vert,
    zmax,
    param_set,
    m;
    zmin = zmin,
    numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
);


t0 = FT(0)
timeend = FT(60)
dt = FT(1)

solver_config =
    ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
mygrid = solver_config.dg.grid;
Q = solver_config.Q;
aux = solver_config.dg.state_auxiliary;
t = ODESolvers.gettime(solver_config.solver)
state_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    Q,
    vars_state_conservative(m, FT),
)
aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    aux,
    vars_state_auxiliary(m, FT);,
)
init_all_vars = OrderedDict(state_vars..., aux_vars...);
init_all_vars["t"] = [t]

ClimateMachine.invoke!(solver_config;);

t = ODESolvers.gettime(solver_config.solver)
state_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    Q,
    vars_state_conservative(m, FT),
)
aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    aux,
    vars_state_auxiliary(m, FT);,
)
final_all_vars = OrderedDict(state_vars..., aux_vars...);
final_all_vars["t"] = [t]


#check that the state at the end matches the state at the beginning within some threshold, everywhere in space.

@test final_all_vars["soil.water.ϑ"] ≈ init_all_vars["soil.water.ϑ"]
=======
m_soil = SoilModel(SoilWaterModel(FT), HeatModel())
sources = ()
m = LandModel(param_set, m_soil, sources)
#test that m.soil == m_soil, just to make sure everything loaded and ran properly.
@test m.soil == m_soil
>>>>>>> master
