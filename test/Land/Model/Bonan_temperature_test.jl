# Test heat equation agrees with solution from Bonan's Supplemental program 5.2
# It doesnt not agree using these parameters, but it does show evolution in time (issue of not changing isnt occuring anymore)
#Not sure re: why this is so slow compared to the heat tutorial. Many other calculations every step, is that it? It's a lot slower.
using MPI
using OrderedCollections
using StaticArrays
using Statistics
using Dierckx
using Test #when finished debugging, can remove. Test is imported in runtests.jl
using Profile

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

#@testset "Bonan temperature test" begin
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
        FT = eltype(state)
        #state.soil.water.ϑ_l = FT(land.soil.water.initialϑ_l(aux))
        #state.soil.water.θ_ice = FT(land.soil.water.initialθ_ice(aux))

        #Density of liquid water (kg/m``^3``)
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

        ϑ_l, θ_ice = get_water_content(land.soil.water, aux, state, time)
        θ_l = volumetric_liquid_fraction(ϑ_l, land.soil.param_functions.porosity)
        c_s = volumetric_heat_capacity(θ_l, θ_ice, land.soil.param_functions.c_ds,
                                       _cp_l, _cp_i)

        state.soil.heat.I = FT(internal_energy(
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
    surface_value = FT(300)
    heat_surface_state = (aux, t) -> surface_value
    # If one wanted a nonzero flux BC, would need to specify entire κ∇T term.
    heat_bottom_flux = (aux, t) -> zero_output
    #This is the initial condition for T. We determine the IC for I using T and θ_ice.
    initial_temp = FT(295.15)
    T_init = (aux) -> initial_temp
    soil_water_model  = PrescribedWaterModel((aux, t) -> zero_output,
                                             (aux, t) -> zero_output
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
    timeend = FT(1)

    # We'll define the time-step based on the [Fourier
    # number](https://en.wikipedia.org/wiki/Fourier_number)
    Δ = min_node_distance(driver_config.grid)

    given_Fourier = FT(0.08);
    Fourier_bound = given_Fourier * Δ^2 / soil_param_functions.κ_dry;
    dt = Fourier_bound

    solver_config =
        ClimateMachine.SolverConfiguration(
            t0,
            timeend,
            driver_config,
            ode_dt = dt
    )
    mygrid = solver_config.dg.grid
    Q = solver_config.Q
    aux = solver_config.dg.state_auxiliary

    @time ClimateMachine.invoke!(solver_config)
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

    #30sec runtime
    # bonan_temperature = reverse([
    #  300.000,
    #  299.961,
    #  299.882,
    #  299.804,
    #  299.725,
    #  299.647,
    #  299.569,
    #  299.491,
    #  299.414,
    #  299.337,
    #  299.261,
    #  299.185,
    #  299.110,
    #  299.035,
    #  298.962,
    #  298.888,
    #  298.816,
    #  298.745,
    #  298.674,
    #  298.604,
    #  298.536,
    #  298.468,
    #  298.402,
    #  298.336,
    #  298.272,
    #  298.209,
    #  298.148,
    #  298.088,
    #  298.029,
    #  297.971,
    #  297.915,
    #  297.861,
    #  297.807,
    #  297.756,
    #  297.706,
    #  297.658,
    #  297.611,
    #  297.566,
    #  297.523,
    #  297.482,
    #  297.442,
    #  297.405,
    #  297.369,
    #  297.334,
    #  297.302,
    #  297.272,
    #  297.244,
    #  297.217,
    #  297.193,
    #  297.170,
    #  297.150,
    #  297.131,
    #  297.115,
    #  297.100,
    #  297.088,
    #  297.078,
    #  297.069,
    #  297.063,
    #  297.059,
    #  297.057
    # ])

    # bonan_z = reverse([
    #    0.000,
    #   -0.847,
    #   -2.542,
    #   -4.237,
    #   -5.932,
    #   -7.627,
    #   -9.322,
    #  -11.017,
    #  -12.712,
    #  -14.407,
    #  -16.102,
    #  -17.797,
    #  -19.492,
    #  -21.186,
    #  -22.881,
    #  -24.576,
    #  -26.271,
    #  -27.966,
    #  -29.661,
    #  -31.356,
    #  -33.051,
    #  -34.746,
    #  -36.441,
    #  -38.136,
    #  -39.831,
    #  -41.525,
    #  -43.220,
    #  -44.915,
    #  -46.610,
    #  -48.305,
    #  -50.000,
    #  -51.695,
    #  -53.390,
    #  -55.085,
    #  -56.780,
    #  -58.475,
    #  -60.169,
    #  -61.864,
    #  -63.559,
    #  -65.254,
    #  -66.949,
    #  -68.644,
    #  -70.339,
    #  -72.034,
    #  -73.729,
    #  -75.424,
    #  -77.119,
    #  -78.814,
    #  -80.508,
    #  -82.203,
    #  -83.898,
    #  -85.593,
    #  -87.288,
    #  -88.983,
    #  -90.678,
    #  -92.373,
    #  -94.068,
    #  -95.763,
    #  -97.458,
    #  -99.153
    # ])

   # #1secruntime
#     bonan_temperature = reverse([
#      300.000,
#      299.768,
#      299.306,
#      298.855,
#      298.419,
#      298.005,
#      297.617,
#      297.259,
#      296.933,
#      296.641,
#      296.382,
#      296.156,
#      295.963,
#      295.798,
#      295.661,
#      295.548,
#      295.456,
#      295.383,
#      295.325,
#      295.280,
#      295.245,
#      295.219,
#      295.199,
#      295.185,
#      295.174,
#      295.167,
#      295.161,
#      295.158,
#      295.155,
#      295.153,
#      295.152,
#      295.151,
#      295.151,
#      295.151,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      295.150,
#      ])

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
#     ])

#     bonan_z = bonan_z ./ 100.0
#     # Create an interpolation from the Bonan data
#     bonan_temperature_continuous = Spline1D(bonan_z, bonan_temperature)
#     bonan_at_clima_z = [bonan_temperature_continuous(i) for i in all_vars["z"]]

    # plot([bonan_at_clima_z all_vars["soil.heat.T"]], all_vars["z"], label = ["Bonan at Clima z" "Clima"])
    # xlabel!("Temperature [K]")
    # ylabel!("Depth [cm]")

    #this is not quite a true L2, because our z values are not equally spaced.
    #MSE = mean((bonan_at_clima_z .- all_vars["soil.heat.T"]) .^ 2.0)
    #@test MSE < 1e-3
#end
