using MPI
using OrderedCollections
using StaticArrays
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
using ClimateMachine
using ClimateMachine.Land
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
ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));

# Load some helper functions for plotting
# include(joinpath(clima_dir, "docs", "plothelpers.jl"));

struct WaterModel end
struct HeatModel end

m_soil = SoilModel(WaterModel(), HeatModel())
sources = ()
m = LandModel(param_set, m_soil, sources)


N_poly = 5;
nelem_vert = 10;
zmax = FT(1);
# driver_config = ClimateMachine.SingleStackConfiguration(
#     "LandModel",
#     N_poly,
#     nelem_vert,
#     zmax,
#     param_set,
#     m,
#     numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
# );

# t0 = FT(0)
# timeend = FT(40)
# dt = FT(0.1)

# solver_config =
#     ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);

# ClimateMachine.invoke!(solver_config);


# function init_state_auxiliary!(land::LandModel, aux::Vars, geom::LocalGeometry)
#     aux.coord = geom.coord
#     init_state_auxiliary!(land, land.soil, aux, geom)
# end


