### Soil water model
export SoilWaterModel
include("./soil_water_parameterizations.jl")
"""
    SoilWaterModel{} <: BalanceLaw

"""

struct SoilWaterModel{FT, IF, VF, MF, HM} <: BalanceLaw
    "Impedance Factor"
    impedance_factor::IF
    "Viscosity Factor"
    viscosity_factor::VF
    "Moisture Factor"
    moisture_factor::MF
    "Hydraulics Model"
    hydraulics::HM # vG, B/C, Haverkamp
end


#Defaults imply a constant hydraulic K = K sat
function SoilWaterModel(::Type{FT};##we need to tell it what FT is in a non-optional way.
                        impedance_factor::AbstractImpedanceFactor{FT} = NoImpedance{FT}(),
                        viscosity_factor::AbstractViscosityFactor{FT} = ConstantViscosity{FT}(),
                        moisture_factor::AbstractMoistureFactor{FT} =  MoistureIndependent{FT}(),
                        hydraulics::AbstractHydraulicsModel{FT} =  vanGenuchten{FT}()
                         
    ) where {FT}
    return SoilWaterModel{FT,##this is a call of the default constructor. To call the outer constructor we dont need to give the type arguments.
                          typeof(impedance_factor),##here we get the concrete types
                          typeof(viscosity_factor),
                          typeof(moisture_factor),
                          typeof(hydraulics)
        }(impedance_factor, viscosity_factor, moisture_factor, hydraulics)
end
