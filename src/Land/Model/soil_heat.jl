#### Soil heat model

abstract type AbstractHeatModel <: BalanceLaw end

function vars_state_conservative(soil::AbstractHeatModel, FT)
    @vars begin end
end

struct ConstantInternalEnergy{FT} <: AbstractHeatModel
    T::FT
end

get_temperature(m::ConstantInternalEnergy) = m.T
