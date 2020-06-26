using MPI
using OrderedCollections
using StaticArrays

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land
using ClimateMachine.Land.SoilWaterParameterizations
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
using Plots
using CSV

FT = Float64;

# When using ClimateMachine, we need to define a function that sets the initial
# state of our model run.

function init_soil_water!(land, state, aux, coordinates, time)
    FT = eltype(state)
    state.soil.water.ν = FT(land.soil.water.initialν)
#    state.ρu = SVector{3, FT}(0, 0, 0) might be a useful ref later for how to initialize vectors.
end;




ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));

 ## Set initial condition

struct HeatModel end

SoilParams = SoilParamSet(porosity = 0.495, Ksat = 0.0443 / (3600*100), S_s = 1e-3)

soil_water_model = SoilWaterModel(FT;
                                  moisture_factor = MoistureDependent{FT}(),
                                  hydraulics = Haverkamp{FT}(),
                                  params = SoilParams,
                                  initialν = 0.24,
                                  surfaceν = 0.494
                                  )
                                  
m_soil = SoilModel(soil_water_model, HeatModel())
sources = ()
m = LandModel(param_set,
              m_soil,
              sources,
              init_soil_water!
              )


N_poly = 5;
nelem_vert = 10;

# Specify the domain boundaries
zmax = FT(0);
zmin = FT(-1)

driver_config = ClimateMachine.SingleStackConfiguration(
    "SoilMoistureModel",
    N_poly,
    nelem_vert,
    zmax,
    param_set,
    m;
    zmin = zmin,
    numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
);

t0 = FT(0)
timeend = FT(60*60*24)

dt = FT(6)

solver_config =
    ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
mygrid = solver_config.dg.grid;
Q = solver_config.Q;
aux = solver_config.dg.state_auxiliary;
all_data = Dict([k => Dict() for k in 0:1]...)
t = ODESolvers.gettime(
    solver_config.solver
)
state_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    Q,
    vars_state_conservative(m, FT),
)
aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    aux,
    vars_state_auxiliary(m, FT);
)
all_vars = OrderedDict(state_vars..., aux_vars...);
all_vars["t"]= [t]
all_data[0] = all_vars
ClimateMachine.invoke!(solver_config);
t = ODESolvers.gettime(
    solver_config.solver
)
state_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    Q,
    vars_state_conservative(m, FT),
)
aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    aux,
    vars_state_auxiliary(m, FT);
)
all_vars = OrderedDict(state_vars..., aux_vars...);
all_vars["t"]= [t]
all_data[1] = all_vars


#compare with haverkamp data here
#need to create interpolation from haverkamp data, then evaluated at our z points
#https://github.com/CliMA/ClimateMachine.jl/blob/agl/ck/soil_water/tutorials/Land/helper_funcs.jl#L154-L190

#wont save the plot in actual test, or even make one.
bonanfile = "./bonan_1day_haverkamp.csv"
bonan_data = CSV.read(bonanfile; header = ["θ","z"],types=Dict(2 => Float64))

plot(all_vars["soil.water.ν"], all_vars["z"]*100, label = "Clima")
plot!(bonan_data[!,1], bonan_data[!,2], label = "Bonan")
savefig("./haverkamp_plot.png")

#bonan_moisture = bonan_data[!,1]
#bonan_z = bonan_data[!,2]


