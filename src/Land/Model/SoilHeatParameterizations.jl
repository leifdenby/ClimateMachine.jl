"""
    SoilHeatParameterizations

Functions for volumetric heat capacity, internal energy as function of temperature,
saturated thermal conductivity, thermal conductivty, the Kersten number, and relative
saturation are included.
"""

module SoilHeatParameterizations

using DocStringExtensions

export AbstractHeatCapacity
    HeatCapacity
    volumetric_heat_capacity
    InternalEnergy
    internal_energy
    ThermalConductivity
    κ_sat
    κ
    KerstenNumber
    Kersten
    S_r

"""
    AbstractHeatCapacity{FT <: AbstractFloat}

Document here.
"""
abstract type AbstractHeatCapacity{FT <: AbstractFloat} end

"""
    HeatCapacity{FT} <: AbstractHeatCapacity{FT}

The necessary parameters for determining heat capacity of soil matrix.
# Fields

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct HeatCapacity{FT} <: AbstractHeatCapacity{FT}
    "Volum. heat capacity liquid water. Units of J m-3 K-1."
    c_l::FT = FT(4.18e6) #get from PlanetParameters
    "Volumetric heat capacity ice. Units of J m-3 K-1."
    c_i::FT = FT(1.93e6) #get from PlanetParameters
end

"""
    volumetric_heat_capacity(
            ϴ_l::FT,
            ϴ_i::FT,
            HeatParams::AbstractHeatCapacity{FT},
            SoilParams::AbstractSoilParameterSet{FT},
            porosity::FT
    ) where {FT}
Compute the expression for volumetric heat capacity.
"""
function volumetric_heat_capacity(
            ϴ_l::FT,
            ϴ_i::FT,
            HeatParams::AbstractHeatCapacity{FT},
            SoilParams::AbstractSoilParameterSet{FT},
            porosity::FT
) where {FT}

    c_l = HeatParams.c_l ## Should the function for ϴ_l be called here or is it ok that it is an input?
    c_i = HeatParams.c_i
    c_ds = SoilParams.c_ds
    c_s = (1 - porosity) * c_ds + ϴ_l * c_l + ϴ_i * c_i
    return c_s
}

"""
    InternalEnergy{FT} <: AbstractInternalEnergy{FT}

The necessary parameters for determining internal energy of
soil matrix as function of temperature.
# Fields

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct InternalEnergy{FT} <: AbstractInternalEnergy{FT}
    "Density of ice. Units of kg m-3."
    ρ_i::FT = FT(916.7) #get from PlanetParameters
    "Arbitrary reference temperature at which the specific internal energy of soil
    and liquid water go to zero. Units of K."
    T_0::FT = FT(273.16) #get from PlanetParameters
    "Specific latent heat of fusion at T_0. Units of J kg-1."
    L_f_0::FT = FT(333.6e3)
end

"""
    internal_energy(
            ϑ_l::FT,
            ϴ_i::FT,
            T::FT,
            params::AbstractInternalEnergy,
    ) where {FT}
Compute the expression for volumetric liquid fraction.
"""
function internal_energy(
            ϑ_l::FT,
            ϴ_i::FT,
            T::FT,
            params::AbstractInternalEnergy{FT}
) where {FT}

    ρ_i = params.ρ_i
    T_0 = params.T_0
    L_f_0 = params.L_f_0
    I = c_s * (T - T0) - ϴ_i * ρ_i * L_f_0
    return
end

"""
    κ_sat(
            ϑ_l::FT,
            ϴ_l::FT,
            ϴ_i::FT,
            porosity::FT,
            SoilParams::AbstractSoilParameterSet{FT}
    ) where {FT}
Compute the expression for saturated thermal conductivity of soil matrix.
"""
function κ_sat(
    θ_l::FT,
    ϴ_l::FT,
    ϴ_i::FT,
    porosity::FT,
    SoilParams::AbstractSoilParameterSet{FT}
) where {FT}

    κ_sat_unfrozen = SoilParams.κ_sat_unfrozen
    κ_sat_frozen = SoilParams.κ_sat_frozen
    ϴ_w = ϴ_l + ϴ_i
    κ_sat = κ_sat_unfrozen^(ϴ_l / ϴ_w) * κ_sat_frozen^(ϴ_i / ϴ_w)
    return κ_sat
end

"""
    κ(
            ϑ_l::FT,
            ϴ_i::FT,
            porosity::FT,
            params::AbstractThermalConductivity
    ) where {FT}
Compute the expression for thermal conductivity of soil matrix.
"""
function κ(
    θ_l::FT,
    ϴ_i::FT,
    porosity::FT,
    SoilParams::AbstractSoilParameterSet{FT}
) where {FT}

    κ_dry = SoilParams.κ_dry
    K_e = Kersten(ϑ, ϴ_i, porosity, params)
    κ_sat = κ_sat(ϑ_l, θ_l, ϴ_i, porosity, c)
    κ = K_e * κ_sat + (1 - K_e) * κ_dry
    return κ
end


"""
    KerstenNumber{FT} <: AbstractKerstenNumber{FT}

The necessary parameters for determining Kersten number.
# Fields

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KerstenNumber{FT} <: AbstractKerstenNumber{FT}
    "Adjustable scale parameter. Unitless."
    a::FT = FT(0.24)
    "Adjustable scale parameter. Unitless."
    b::FT = FT(18.1)
end
"""
    Kersten(
            ϑ_l::FT,
            ϴ_i::FT,
            porosity::FT,
            params::AbstractKerstenNumber{FT}
    ) where {FT}
Compute the expression for the Kersten number.
"""
function Kersten(
    θ_l::FT,
    ϴ_i::FT,
    porosity::FT,
    HeatParams::AbstractHeatCapacity{FT},
    SoilParams::AbstractSoilParameterSet{FT}
) where {FT}

    a = HeatParams.a
    b = HeatParams.a
    ν_om = SoilParams.ν_om
    ν_sand = SoilParams.ν_sand
    ν_gravel = SoilParams.ν_gravel
    S_r = S_r(ϑ_l, θ_i, porosity)
    if ϴ_i = 0
        K_e = S_r^((1 + ν_om - a * ν_sand - ν_gravel) / 2)*([1 + exp(-b * S_r)]^(-3) - ((1 - S_r) / 2)^3)^(1 - ν_om)
    else
        K_e=S_r^(1 + ν_om)
    end
    return K_e
end

"""
    S_r(
            ϑ_l::FT,
            ϴ_i::FT,
            porosity::FT
    ) where {FT}
Compute the expression for relative saturation.
"""
function S_r(
    θ_l::FT,
    ϴ_i::FT,
    porosity::FT
) where {FT}

    S_r=(ϑ_l + ϴ_i) / porosity
    return S_r
end

end # Module
