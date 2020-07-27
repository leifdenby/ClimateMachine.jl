# Test heat equation agrees with solution from Bonan's Supplemental program 5.2
using MPI
using OrderedCollections
using StaticArrays
using Statistics
using Dierckx

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land
using ClimateMachine.Land.SoilWaterParameterizations
using ClimateMachine.Land.SoilHeatParameterizations
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

# This has the interpolation functions if we went with Interpolations.jl
#include("./helperfunc.jl")

FT = Float64;

# Density of liquid water (kg/m``^3``)
_ρ_l = FT(ρ_cloud_liq(param_set))
# Density of ice water (kg/m``^3``)
_ρ_i = FT(ρ_cloud_ice(param_set))
# Volum. isoboric heat capacity liquid water (J/m3/K)
_cp_l = FT(cp_l(param_set) * _ρ_l)
# Volumetric isoboric heat capacity ice (J/m3/K)
_cp_i = FT(cp_i(param_set) * _ρ_i)
# Density of ice water (kg/m``^3``)
_ρ_i = FT(ρ_cloud_ice(param_set))
# Reference temperature (K)
_T_ref = FT(T_0(param_set))
# Latent heat of fusion at ``T_0`` (J/kg)
_LH_f0 = FT(LH_f0(param_set))

# When using ClimateMachine, we need to define a function that sets the initial
# state of our model run.

function init_soil!(land, state, aux, coordinates, time)
    FT = eltype(state)
    #state.soil.water.ϑ_l = FT(land.soil.water.initialϑ_l(aux))
    #state.soil.water.θ_ice = FT(land.soil.water.initialθ_ice(aux))
    state.soil.heat.I = FT(land.soil.heat.params.ρc * land.soil.heat.initialT(aux)) # land.soil.heat.initialT(aux))
    #    state.ρu = SVector{3, FT}(0, 0, 0) might be a useful ref later for how to initialize vectors.
end;

ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));

SoilParams = SoilParamSet(
        porosity = FT(0.495),
        Ksat = 0.0443 / (3600*100),
        S_s = 1e-3,
        ν_gravel = 0.1,
        ν_om = 0.1,
        ν_sand = 0.1,
        c_ds = 1e6,
        κ_dry = 1.5,
        κ_sat_unfrozen = 0.57,
        κ_sat_frozen = 2.29,
        ρc = 1.0,
        α = 0.01,
        a = 0.24,
        b = 18.1
        )
# Keep in mind that what is passed is aux⁻
# Fluxes are multiplied by ẑ (normal to the surface, -normal to the bottom,
# where normal points out of the domain.)

# water_surface_state = (aux, t) -> FT(0.494)
# water_bottom_flux = (aux, t) -> FT(aux.soil.water.κ * 1.0)
# ϑ_l0 = (aux) -> FT(0.24)

heat_surface_state = (aux, t) -> FT(300)
heat_bottom_flux = (aux, t) -> FT(0)
T_init = (aux) -> FT(295.15)

#soil_water_model = PrescribedWaterModel{FT}()

soil_water_model = PrescribedWaterModel(
    FT;
    ϑ_l = FT(0.0),
    θ_ice = FT(0.0)
    )

# soil_water_model = SoilWaterModel(
#     FT;
#     moisture_factor = MoistureDependent{FT}(),
#     hydraulics = Haverkamp{FT}(),
#     params = SoilParams,
#     initialϑ_l = ϑ_l0,
#     dirichlet_bc = Dirichlet(
#         surface_state = water_surface_state,
#         bottom_state = nothing,
#     ),
#     neumann_bc = Neumann(
#         surface_flux = nothing,
#         bottom_flux = water_bottom_flux
#     ),
# )

#soil_heat_model = PrescribedTemperatureModel{FT}()
soil_heat_model = SoilHeatModel(
    FT;
    params = SoilParams,
    initialT = T_init,
    dirichlet_bc = Dirichlet(
        surface_state = heat_surface_state,
        bottom_state = nothing
    ),
    neumann_bc = Neumann(
        surface_flux = nothing,
        bottom_flux = heat_bottom_flux
    ),
)

m_soil = SoilModel(soil_water_model, soil_heat_model)
sources = ()
m = LandModel(
    param_set,
    m_soil;
    source = sources,
    init_state_conservative = init_soil!,
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

# We'll define the time-step based on the [Fourier
# number](https://en.wikipedia.org/wiki/Fourier_number)
Δ = min_node_distance(driver_config.grid)

given_Fourier = FT(0.08);
Fourier_bound = given_Fourier * Δ^2 / SoilParams.α;
dt = Fourier_bound

solver_config =
    ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
mygrid = solver_config.dg.grid;
Q = solver_config.Q;
aux = solver_config.dg.state_auxiliary;

ClimateMachine.invoke!(solver_config);
t = ODESolvers.gettime(solver_config.solver)
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
all_vars["t"] = [t]
#plot(all_vars["soil.heat.T"],all_vars["z"])

bonan_temperature = reverse([
     300.000,
     299.989,
     299.966,
     299.944,
     299.921,
     299.899,
     299.877,
     299.854,
     299.832,
     299.810,
     299.787,
     299.765,
     299.743,
     299.721,
     299.699,
     299.677,
     299.655,
     299.633,
     299.612,
     299.590,
     299.569,
     299.547,
     299.526,
     299.505,
     299.484,
     299.463,
     299.442,
     299.422,
     299.401,
     299.381,
     299.361,
     299.341,
     299.321,
     299.301,
     299.281,
     299.262,
     299.243,
     299.224,
     299.205,
     299.187,
     299.168,
     299.150,
     299.132,
     299.114,
     299.097,
     299.079,
     299.062,
     299.046,
     299.029,
     299.012,
     298.996,
     298.980,
     298.965,
     298.949,
     298.934,
     298.919,
     298.905,
     298.891,
     298.877,
     298.863,
     298.849,
     298.836,
     298.823,
     298.810,
     298.798,
     298.786,
     298.774,
     298.763,
     298.752,
     298.741,
     298.730,
     298.720,
     298.710,
     298.701,
     298.692,
     298.683,
     298.674,
     298.666,
     298.658,
     298.650,
     298.643,
     298.636,
     298.629,
     298.623,
     298.617,
     298.612,
     298.606,
     298.601,
     298.597,
     298.593,
     298.589,
     298.585,
     298.582,
     298.579,
     298.577,
     298.575,
     298.573,
     298.572,
     298.570,
     298.570,
     298.569,
])

bonan_z = reverse([
       0.000,
      -0.500,
      -1.500,
      -2.500,
      -3.500,
      -4.500,
      -5.500,
      -6.500,
      -7.500,
      -8.500,
      -9.500,
     -10.500,
     -11.500,
     -12.500,
     -13.500,
     -14.500,
     -15.500,
     -16.500,
     -17.500,
     -18.500,
     -19.500,
     -20.500,
     -21.500,
     -22.500,
     -23.500,
     -24.500,
     -25.500,
     -26.500,
     -27.500,
     -28.500,
     -29.500,
     -30.500,
     -31.500,
     -32.500,
     -33.500,
     -34.500,
     -35.500,
     -36.500,
     -37.500,
     -38.500,
     -39.500,
     -40.500,
     -41.500,
     -42.500,
     -43.500,
     -44.500,
     -45.500,
     -46.500,
     -47.500,
     -48.500,
     -49.500,
     -50.500,
     -51.500,
     -52.500,
     -53.500,
     -54.500,
     -55.500,
     -56.500,
     -57.500,
     -58.500,
     -59.500,
     -60.500,
     -61.500,
     -62.500,
     -63.500,
     -64.500,
     -65.500,
     -66.500,
     -67.500,
     -68.500,
     -69.500,
     -70.500,
     -71.500,
     -72.500,
     -73.500,
     -74.500,
     -75.500,
     -76.500,
     -77.500,
     -78.500,
     -79.500,
     -80.500,
     -81.500,
     -82.500,
     -83.500,
     -84.500,
     -85.500,
     -86.500,
     -87.500,
     -88.500,
     -89.500,
     -90.500,
     -91.500,
     -92.500,
     -93.500,
     -94.500,
     -95.500,
     -96.500,
     -97.500,
     -98.500,
     -99.500,
])

bonan_z = bonan_z ./ 100.0
# Create an interpolation from the Bonan data
bonan_temperature_continuous = Spline1D(bonan_z, bonan_temperature)
bonan_at_clima_z = [bonan_temperature_continuous(i) for i in all_vars["z"]]

# plot([bonan_at_clima_z all_vars["soil.heat.T"]], all_vars["z"], label = ["Bonan at Clima z" "Clima"])
# xlabel!("Temperature [K]")
# ylabel!("Depth [cm]")

#this is not quite a true L2, because our z values are not equally spaced.
MSE = mean((bonan_at_clima_z .- all_vars["soil.heat.T"]) .^ 2.0)
#@test MSE < 1e-5
