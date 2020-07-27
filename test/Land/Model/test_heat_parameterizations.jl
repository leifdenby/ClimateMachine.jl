using MPI
using OrderedCollections
using StaticArrays
using Test

using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land.SoilHeatParameterizations

FT = Float64;
# Density of liquid water (kg/m``^3``)
_ρ_l = FT(ρ_cloud_liq(param_set))
# Density of ice water (kg/m``^3``)
_ρ_i = FT(ρ_cloud_ice(param_set))
# Volum. isoboric heat capacity liquid water (J/m3/K)
_cp_l = FT(cp_l(param_set) * _ρ_l)
# Volumetric isoboric heat capacity ice (J/m3/K)
_cp_i = FT(cp_i(param_set) * _ρ_i)
# Reference temperature (K)
_T_ref = FT(T_0(param_set))
# Latent heat of fusion at ``T_0`` (J/kg)
_LH_f0 = FT(LH_f0(param_set))

# Internal tests
# top & bc bcs equal same constant

# Unit tests
# Haverkamp just heat
# E conservation ## for future & freeze thaw (two separate branches)
# Compare against analytical soln (see soln on overleaf)
# Other unit test for both water and heat

@testset "Land heat parameterizations" begin
    @test volumetric_heat_capacity(
        0.25,
        0.05,
        1e6,
        _cp_l,
        _cp_i
    ) == FT(2.1415e+06)

    @test internal_energy(
        0.05,
        2.1415e+06,
        300,
        _T_ref,
        _ρ_i,
        _LH_f_0
    ) == FT(4.2187e+07)

    @test Saturated_thermal_conductivity(
        0.25,
        0.05,
        0.57,
        2.29
    ) == FT(0.7187)

    @test Relative_saturation(
        0.25,
        0.05,
        0.4
    ) == FT(0.75)

    # Test branching in Kersten_Number

    # ice fraction = 0
    @test Kersten_Number(
        0.0,
        0.75,
        0.24,
        18.1,
        0.1,
        0.1,
        0.1
    ) == FT(0.8675)

    # ice fraction ~= 0
    @test Kersten_Number(
        0.05,
        0.75,
        0.24,
        18.1,
        0.1,
        0.1,
        0.1
    ) == FT(0.7287)

    @test Thermal_conductivity(
        1.5,
        0.7287,
        0.7187
    ) == FT(0.8222)

    # @test_throws DomainError effective_saturation(0.5, -1.0)

    # to test an array example
    # test_array = [0.5, 1.0]
    # @test Kersten_Number(test_array, Ref(S_r), Ref(a), Ref(b), Ref(ν_om), Ref(ν_sand), Ref(ν_gravel)) ≈
    #       .S_r.^((1 + ν_om - a * ν_sand - ν_gravel) / 2).*([1 + exp(-b .* S_r)].^(-3) .- ((1 .- S_r) ./ 2).^3).^(1 - ν_om)
end
