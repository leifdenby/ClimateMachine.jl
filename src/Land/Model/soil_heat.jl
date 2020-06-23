### Soil heat model
export SoilHeatModel

"""
    SoilHeatModel <: BalanceLaw

The balance law for internal energy in soil.

"""

abstract type AbstractHeatModel <: BalanceLaw end

struct SoilHeatModel <: AbstractHeatModel end

"""
    vars_state_conservative(soil::AbstractHeatModel, FT)

"""

function vars_state_conservative(soil::AbstractHeatModel, FT)
    @vars begin end
end

struct ConstantInternalEnergy{FT} <: AbstractHeatModel
    T::FT
end

"""
    get_temperature(m::ConstantInternalEnergy)

"""

get_temperature(m::ConstantInternalEnergy) = m.T
