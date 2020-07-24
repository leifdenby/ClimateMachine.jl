"""
    SoilHeatParameterizations

Functions for volumetric heat capacity, internal energy as function of temperature,
saturated thermal conductivity, thermal conductivty, the Kersten number, and relative
saturation are included.
"""

module SoilHeatParameterizations

using DocStringExtensions

using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, cp_l, cp_i, T_0, LH_f0
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

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

export volumetric_heat_capacity,
    internal_energy,
    Saturated_thermal_conductivity,
    Thermal_conductivity,
    Relative_saturation,
    Kersten_Number

"""
    volumetric_heat_capacity(
        ϴ_l::FT,
        ϴ_i::FT,
        porosity::FT,
        c_ds::FT,
        cp_l::FT,
        cp_i::FT
    ) where {FT}
Compute the expression for volumetric heat capacity.
"""
function volumetric_heat_capacity(
    ϴ_l::FT,
    ϴ_i::FT,
    porosity::FT,
    c_ds::FT,
    cp_l::FT,
    cp_i::FT
) where {FT}

    c_s = (1 - porosity) * c_ds + ϴ_l *cp_l + ϴ_i * cp_i
    return c_s
end

"""
    internal_energy(
        ϑ_l::FT,
        ϴ_i::FT,
        c_s::FT,
        T::FT
        T_ref::FT,
        ρ_i::FT,
        LH_f_0::FT
    ) where {FT}
Compute the expression for volumetric liquid fraction.
"""
function internal_energy(
    ϑ_l::FT,
    ϴ_i::FT,
    c_s::FT,
    T::FT,
    T_ref::FT,
    ρ_i::FT,
    LH_f_0::FT
) where {FT}

    I = c_s * (T - T_ref) - ϴ_i * ρ_i * LH_f_0
    return
end

"""
    Saturated_thermal_conductivity(
        ϑ_l::FT,
        ϴ_l::FT,
        ϴ_i::FT,
        porosity::FT,
        κ_sat_unfrozen::FT,
        κ_sat_frozen::FT
    ) where {FT}
Compute the expression for saturated thermal conductivity of soil matrix.
"""
function Saturated_thermal_conductivity(
    θ_l::FT,
    ϴ_l::FT,
    ϴ_i::FT,
    porosity::FT,
    κ_sat_unfrozen::FT,
    κ_sat_frozen::FT
) where {FT}

    ϴ_w = ϴ_l + ϴ_i
    κ_sat = κ_sat_unfrozen^(ϴ_l / ϴ_w) * κ_sat_frozen^(ϴ_i / ϴ_w)
    return κ_sat
end

"""
    Thermal_conductivity(
        ϑ_l::FT,
        θ_l::FT,
        ϴ_i::FT,
        porosity::FT,
        κ_dry::FT,
        K_e::FT,
        κ_sat::FT
    ) where {FT}
Compute the expression for thermal conductivity of soil matrix.
"""
function Thermal_conductivity(
    ϑ_l::FT,
    ϴ_l::FT,
    ϴ_i::FT,
    porosity::FT,
    κ_dry::FT,
    K_e::FT,
    κ_sat::FT
) where {FT}

    κ = K_e * κ_sat + (1 - K_e) * κ_dry
    return κ
end

"""
    Relative_saturation(
            ϑ_l::FT,
            ϴ_i::FT,
            porosity::FT
    ) where {FT}
Compute the expression for relative saturation.
"""
function Relative_saturation(
    θ_l::FT,
    ϴ_i::FT,
    porosity::FT
) where {FT}

    S_r=(ϑ_l + ϴ_i) / porosity
    return S_r
end

"""
    Kersten_Number(
        θ_l::FT,
        ϴ_i::FT,
        porosity::FT,
        S_r::FT,
        a::FT,
        b::FT,
        ν_om::FT,
        ν_sand::FT,
        ν_gravel::FT
    ) where {FT}
Compute the expression for the Kersten number.
"""
function Kersten_Number(
    θ_l::FT,
    ϴ_i::FT,
    porosity::FT,
    S_r::FT,
    a::FT,
    b::FT,
    ν_om::FT,
    ν_sand::FT,
    ν_gravel::FT
) where {FT}

    if ϴ_i == 0 # This might give an error due to it not being exactly equal to 0?
        K_e = S_r^((1 + ν_om - a * ν_sand - ν_gravel) / 2)*([1 + exp(-b * S_r)]^(-3) - ((1 - S_r) / 2)^3)^(1 - ν_om)
    else
        K_e = S_r^(1 + ν_om)
    end
    return K_e
end

end # Module
