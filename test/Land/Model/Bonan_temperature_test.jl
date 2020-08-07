# Test heat equation agrees with solution from Bonan's Supplemental program 5.2
# It doesnt not agree using these parameters, but it does show evolution in time (issue of not changing isnt occuring anymore)
#Not sure re: why this is so slow compared to the heat tutorial. Many other calculations every step, is that it? It's a lot slower.
using MPI
using OrderedCollections
using StaticArrays
using Statistics
using Dierckx
using Test

using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, cp_l, cp_i, T_0, LH_f0
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
using ClimateMachine.BalanceLaws:
    BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux, vars_state
import ClimateMachine.DGMethods: calculate_dt

# This has the interpolation functions if we went with Interpolations.jl
#include("./helperfunc.jl")

function calculate_dt(dg, model::LandModel, Q, Courant_number, t, direction)
    Δt = one(eltype(Q))
    CFL = DGMethods.courant(diffusive_courant, dg, model, Q, Δt, t, direction)
    return Courant_number / CFL
end
function diffusive_courant(
    m::LandModel,
    state::Vars,
    aux::Vars,
    diffusive::Vars,
    Δx,
    Δt,
    t,
    direction,
)
    return Δt * m.soil.param_functions.κ_dry / (Δx * Δx)
end

@testset "Bonan temperature test" begin
    ClimateMachine.init()
    mpicomm = MPI.COMM_WORLD

    FT = Float32

    # # Density of liquid water (kg/m``^3``)
    # _ρ_l = FT(ρ_cloud_liq(param_set))
    # # Density of ice water (kg/m``^3``)
    # _ρ_i = FT(ρ_cloud_ice(param_set))
    # # Volum. isoboric heat capacity liquid water (J/m3/K)
    # _cp_l = FT(cp_l(param_set) * _ρ_l)
    # # Volumetric isoboric heat capacity ice (J/m3/K)
    # _cp_i = FT(cp_i(param_set) * _ρ_i)
    # # Density of ice water (kg/m``^3``)
    # _ρ_i = FT(ρ_cloud_ice(param_set))
    # # Reference temperature (K)
    # _T_ref = FT(T_0(param_set))
    # # Latent heat of fusion at ``T_0`` (J/kg)
    # _LH_f0 = FT(LH_f0(param_set))


    function init_soil!(land, state, aux, coordinates, time)
        myfloat = eltype(state)
        #state.soil.water.ϑ_l = FT(land.soil.water.initialϑ_l(aux))
        #state.soil.water.θ_ice = FT(land.soil.water.initialθ_ice(aux))

        #Density of liquid water (kg/m``^3``)
        _ρ_l = myfloat(ρ_cloud_liq(param_set))
        # Density of ice water (kg/m``^3``)
        _ρ_i = myfloat(ρ_cloud_ice(param_set))
        # Volum. isoboric heat capacity liquid water (J/m3/K)
        _cp_l = myfloat(cp_l(param_set) * _ρ_l)
        # Volumetric isoboric heat capacity ice (J/m3/K)
        _cp_i = myfloat(cp_i(param_set) * _ρ_i)
        # Density of ice water (kg/m``^3``)
        _ρ_i = myfloat(ρ_cloud_ice(param_set))
        # Reference temperature (K)
        _T_ref = myfloat(T_0(param_set))
        # Latent heat of fusion at ``T_0`` (J/kg)
        _LH_f0 = myfloat(LH_f0(param_set))

       ϑ_l, θ_ice = get_water_content(land.soil.water, aux, state, time)
       θ_l = volumetric_liquid_fraction(ϑ_l, land.soil.param_functions.porosity)
       c_s = volumetric_heat_capacity(θ_l, θ_ice, land.soil.param_functions.c_ds,
                                      _cp_l, _cp_i)

        state.soil.heat.I = myfloat(internal_energy(
        θ_ice,
        c_s,
        land.soil.heat.initialT(aux),
        _T_ref,
        _ρ_i,
        _LH_f0))


    end

    soil_param_functions = SoilParamFunctions{FT}(
            porosity = 0.495,
            Ksat = 0.0443 / (3600*100),
            S_s = 1e-3,
            ν_gravel = 0.1,
            ν_om = 0.1,
            ν_sand = 0.1,
            c_ds = 1,#e6,
            κ_dry = 0.01,#1.5,
            κ_sat_unfrozen = 0.57,
            κ_sat_frozen = 2.29,
        #α/ρc = κ. To run previous test with α = 0.01 and ρc = 1, choose κ_dry = 0.01, and cds = 1, and make sure
        #the water content is 0 for ice and liquid for all time/space. ρc and α are no longer used.
            a = 0.24,
            b = 18.1
            )
    # Keep in mind that what is passed is aux⁻
    # Fluxes are multiplied by ẑ (normal to the surface, -normal to the bottom,
    # where normal points out of the domain.)

    # These are boundary conditions on ϑ_l and/or K∇h. We do not need boundary
    # conditions for θ_ice, since we only solve an ODE for it.
    # water_surface_state = (aux, t) -> FT(0.494)
    # water_bottom_flux = (aux, t) -> FT(aux.soil.water.K * 1.0)

    #Initial condition for ϑ. The default for ice is 0 for all time/space.
    # ϑ_l0 = (aux) -> FT(0.24)

    #Specify boundary condition on T and/or on κ∇T.
    #the BC on T is converted to a BC on I inside the source code.
    zero_output = FT(0.0)
    bottom_value = FT(300)
    heat_surface_flux = (aux, t) -> zero_output
    # If one wanted a nonzero flux BC, would need to specify entire κ∇T term.
    heat_bottom_state = (aux, t) -> bottom_value
    #This is the initial condition for T. We determine the IC for I using T and θ_ice.
    initial_temp = FT(295.15)
    T_init = (aux) -> initial_temp

    soil_water_model  = PrescribedWaterModel((aux, t) -> zero_output,
                                             (aux, t) -> zero_output
                                             )

    soil_heat_model = SoilHeatModel(
        FT;
        initialT = T_init,
        dirichlet_bc = Dirichlet(
            surface_state = nothing,
            bottom_state = heat_bottom_state
        ),
        neumann_bc = Neumann(
            surface_flux = heat_surface_flux,
            bottom_flux = nothing
        ),
    )

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    sources = ()
    m = LandModel(
        param_set,
        m_soil;
        source = sources,
        init_state_prognostic = init_soil!,
    )

    N_poly = 5
    nelem_vert = 10

    # Specify the domain boundaries
    zmax = FT(0)
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
    )

    t0 = FT(0)
    timeend = FT(40)

    solver_config =
        ClimateMachine.SolverConfiguration(
            t0,
            timeend,
            driver_config;
            Courant_number = FT(0.7),
            CFL_direction = VerticalDirection(),
    )
    mygrid = solver_config.dg.grid
    Q = solver_config.Q
    aux = solver_config.dg.state_auxiliary

    ClimateMachine.invoke!(solver_config)
    t = ODESolvers.gettime(solver_config.solver)
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        Q,
        vars_state(m, Prognostic(), FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        aux,
        vars_state(m, Auxiliary(), FT),
    )
    grad_vars = SingleStackUtils.get_vars_from_nodal_stack(
    mygrid,
    solver_config.dg.state_gradient_flux,
    vars_state(m, GradientFlux(), FT),
    )
    all_vars = OrderedDict(state_vars..., aux_vars...)
    all_vars["t"] = [t]
    #plot(all_vars["soil.heat.T"],all_vars["z"])

    # bonan_temperature = reverse([
    #      300.000,
    #      299.989,
    #      299.966,
    #      299.944,
    #      299.921,
    #      299.899,
    #      299.877,
    #      299.854,
    #      299.832,
    #      299.810,
    #      299.787,
    #      299.765,
    #      299.743,
    #      299.721,
    #      299.699,
    #      299.677,
    #      299.655,
    #      299.633,
    #      299.612,
    #      299.590,
    #      299.569,
    #      299.547,
    #      299.526,
    #      299.505,
    #      299.484,
    #      299.463,
    #      299.442,
    #      299.422,
    #      299.401,
    #      299.381,
    #      299.361,
    #      299.341,
    #      299.321,
    #      299.301,
    #      299.281,
    #      299.262,
    #      299.243,
    #      299.224,
    #      299.205,
    #      299.187,
    #      299.168,
    #      299.150,
    #      299.132,
    #      299.114,
    #      299.097,
    #      299.079,
    #      299.062,
    #      299.046,
    #      299.029,
    #      299.012,
    #      298.996,
    #      298.980,
    #      298.965,
    #      298.949,
    #      298.934,
    #      298.919,
    #      298.905,
    #      298.891,
    #      298.877,
    #      298.863,
    #      298.849,
    #      298.836,
    #      298.823,
    #      298.810,
    #      298.798,
    #      298.786,
    #      298.774,
    #      298.763,
    #      298.752,
    #      298.741,
    #      298.730,
    #      298.720,
    #      298.710,
    #      298.701,
    #      298.692,
    #      298.683,
    #      298.674,
    #      298.666,
    #      298.658,
    #      298.650,
    #      298.643,
    #      298.636,
    #      298.629,
    #      298.623,
    #      298.617,
    #      298.612,
    #      298.606,
    #      298.601,
    #      298.597,
    #      298.593,
    #      298.589,
    #      298.585,
    #      298.582,
    #      298.579,
    #      298.577,
    #      298.575,
    #      298.573,
    #      298.572,
    #      298.570,
    #      298.570,
    #      298.569,
    # ])
    
    # bonan_z = reverse([
    #        0.000,
    #       -0.500,
    #       -1.500,
    #       -2.500,
    #       -3.500,
    #       -4.500,
    #       -5.500,
    #       -6.500,
    #       -7.500,
    #       -8.500,
    #       -9.500,
    #      -10.500,
    #      -11.500,
    #      -12.500,
    #      -13.500,
    #      -14.500,
    #      -15.500,
    #      -16.500,
    #      -17.500,
    #      -18.500,
    #      -19.500,
    #      -20.500,
    #      -21.500,
    #      -22.500,
    #      -23.500,
    #      -24.500,
    #      -25.500,
    #      -26.500,
    #      -27.500,
    #      -28.500,
    #      -29.500,
    #      -30.500,
    #      -31.500,
    #      -32.500,
    #      -33.500,
    #      -34.500,
    #      -35.500,
    #      -36.500,
    #      -37.500,
    #      -38.500,
    #      -39.500,
    #      -40.500,
    #      -41.500,
    #      -42.500,
    #      -43.500,
    #      -44.500,
    #      -45.500,
    #      -46.500,
    #      -47.500,
    #      -48.500,
    #      -49.500,
    #      -50.500,
    #      -51.500,
    #      -52.500,
    #      -53.500,
    #      -54.500,
    #      -55.500,
    #      -56.500,
    #      -57.500,
    #      -58.500,
    #      -59.500,
    #      -60.500,
    #      -61.500,
    #      -62.500,
    #      -63.500,
    #      -64.500,
    #      -65.500,
    #      -66.500,
    #      -67.500,
    #      -68.500,
    #      -69.500,
    #      -70.500,
    #      -71.500,
    #      -72.500,
    #      -73.500,
    #      -74.500,
    #      -75.500,
    #      -76.500,
    #      -77.500,
    #      -78.500,
    #      -79.500,
    #      -80.500,
    #      -81.500,
    #      -82.500,
    #      -83.500,
    #      -84.500,
    #      -85.500,
    #      -86.500,
    #      -87.500,
    #      -88.500,
    #      -89.500,
    #      -90.500,
    #      -91.500,
    #      -92.500,
    #      -93.500,
    #      -94.500,
    #      -95.500,
    #      -96.500,
    #      -97.500,
    #      -98.500,
    #      -99.500,
    # ])
bonan_z = reverse([
       0.000,
      -0.847,
      -2.542,
      -4.237,
      -5.932,
      -7.627,
      -9.322,
     -11.017,
     -12.712,
     -14.407,
     -16.102,
     -17.797,
     -19.492,
     -21.186,
     -22.881,
     -24.576,
     -26.271,
     -27.966,
     -29.661,
     -31.356,
     -33.051,
     -34.746,
     -36.441,
     -38.136,
     -39.831,
     -41.525,
     -43.220,
     -44.915,
     -46.610,
     -48.305,
     -50.000,
     -51.695,
     -53.390,
     -55.085,
     -56.780,
     -58.475,
     -60.169,
     -61.864,
     -63.559,
     -65.254,
     -66.949,
     -68.644,
     -70.339,
     -72.034,
     -73.729,
     -75.424,
     -77.119,
     -78.814,
     -80.508,
     -82.203,
     -83.898,
     -85.593,
     -87.288,
     -88.983,
     -90.678,
     -92.373,
     -94.068,
     -95.763,
     -97.458,
     -99.153
    ])

    bonan_temperature = reverse([
     300.000,
     299.866,
     299.599,
     299.333,
     299.071,
     298.813,
     298.561,
     298.315,
     298.077,
     297.847,
     297.626,
     297.416,
     297.215,
     297.026,
     296.847,
     296.680,
     296.523,
     296.378,
     296.244,
     296.121,
     296.009,
     295.906,
     295.813,
     295.729,
     295.654,
     295.586,
     295.526,
     295.473,
     295.427,
     295.386,
     295.350,
     295.319,
     295.292,
     295.269,
     295.249,
     295.233,
     295.218,
     295.206,
     295.196,
     295.188,
     295.181,
     295.175,
     295.170,
     295.166,
     295.163,
     295.160,
     295.158,
     295.156,
     295.155,
     295.154,
     295.153,
     295.152,
     295.152,
     295.151,
     295.151,
     295.151,
     295.151,
     295.151,
     295.150,
     295.150
     ])

#    bonan_z = bonan_z ./ 100.0
    # Create an interpolation from the Bonan data
#    bonan_temperature_continuous = Spline1D(bonan_z, bonan_temperature)
#    bonan_at_clima_z = [bonan_temperature_continuous(i) for i in all_vars["z"]]

    # plot([bonan_at_clima_z all_vars["soil.heat.T"]], all_vars["z"], label = ["Bonan at Clima z" "Clima"])
    # xlabel!("Temperature [K]")
    # ylabel!("Depth [cm]")

    #this is not quite a true L2, because our z values are not equally spaced.
#    MSE = mean((bonan_at_clima_z .- all_vars["soil.heat.T"]) .^ 2.0)
    #@test MSE < 1e-3
end
