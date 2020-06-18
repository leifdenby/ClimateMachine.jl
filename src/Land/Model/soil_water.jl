#### Soil water model

export SoilWaterModel

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


abstract type AbstractImpedanceFactor{FT<:AbstractFloat} end
abstract type AbstractViscosityFactor{FT<:AbstractFloat} end
abstract type AbstractMoistureFactor{FT<:AbstractFloat} end
"""
    Hydraulics model is used in the moisture factor in hydraulic conductivity and in the matric potential. The single hydraulics model choice sets both of these.
"""
abstract type AbstractHydraulicsModel{FT <: AbstractFloat} end

"""
    vanGenuchten{FT} <: AbstractHydraulicsModel{FT}

Temporary stand-in for a hydraulics model. Test case for debugging before adding in full vG, Haverkamp, and Brooks and Corey.
     
"""

Base.@kwdef struct vanGenuchten{FT} <: AbstractHydraulicsModel{FT}
    fake_constant::FT = FT(5.0)
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



struct MoistureIndependent{FT} <: AbstractMoistureFactor{FT}
end

struct MoistureDependent{FT} <: AbstractMoistureFactor{FT}
end


function moisture_factor(
    mm::MoistureDependent{FT},
    hm::vanGenuchten{FT},
    S_l::FT
)where {FT}
    coefficient = hm.fake_constant
    my_standin_function = FT(S_l*coefficient)
    return my_standin_function
end


    
function moisture_factor(
    mm::MoistureIndependent{FT}
) where {FT}
    Factor = FT(1.0)
    return Factor
end


struct ConstantViscosity{FT} <: AbstractViscosityFactor{FT}
end

Base.@kwdef struct TemperatureDependentViscosity{FT} <: AbstractViscosityFactor{FT}
    γ::FT = FT(2.64e-2)
    T_ref::FT = FT(288.0)   
end

function viscosity_factor(
    vm::ConstantViscosity{FT}
) where {FT}
    Theta = FT(1.0)
    return Theta
end

#function viscosity_factor(
#    vm::TemperatureDependentViscosity{FT},
#    T::FT
#) where {FT}
#    γ = vm.γ
#    T_ref = vm.T_ref
#    factor = FT(γ*(T-T_ref))
#    Theta = FT(exp(factor))
#    return Theta
#end

struct NoImpedance{FT} <: AbstractImpedanceFactor{FT}
end

Base.@kwdef struct IceImpedance{FT} <: AbstractImpedanceFactor{FT}
    Ω::FT = FT(7)
end

function impedance_factor(
    imp::NoImpedance{FT}
) where {FT}
    gamma = FT(1.0)
    return gamma
end

#function impedance_factor(
#    imp::IceImpedance{FT},
#    θ_ice::FT,
#    porosity::FT
#) where {FT}
#    Ω = imp.Ω
#    S_ice = θ_ice/porosity 
#    gamma = FT(10.0^(-Ω*S_ice))
#    return gamma
#end

function hydraulic_conductivity(
    impedance::NoImpedance{FT},
    viscosity::ConstantViscosity{FT},
    moisture::MoistureIndependent{FT}
) where {FT}
    K = FT(viscosity_factor(viscosity)*impedance_factor(impedance)*moisture_factor())
    return K
end
    


#function hydraulic_conductivity{FT}(cm::ZeroConductivity) where {FT}
#    return 0
#end
    

#struct ZeroConductivity{FT} <: AbstractConductivityModel{FT}
#end

#struct ConductiveModel{FT, IF, VF, MF} <: AbstractConductivityModel{FT}##these are all concrete types
#    impedance_factor::IF
#    viscosity_factor::VF
#    moisture_factor::MF
#end

#Defaults imply a constant hydraulic K = K sat
#function ConductiveModel(::Type{FT};##we need to tell it what FT is in a non-optional way.
#                         impedance_factor::AbstractImpedanceFactor{FT} = NoImpedance{FT}(),
#                         viscosity_factor::AbstractViscosityFactor{FT} = ConstantViscosity{FT}(),
#                         moisture_factor::AbstractMoistureFactor{FT} =  MoistureIndependent{FT}()
#                         
#    ) where {FT}
#    return ConductiveModel{FT,##this is a call of the default constructor. To call the outer constructor we dont need to give the type arguments.
#        typeof(impedance_factor),##here we get the concrete types
#        typeof(viscosity_factor),
#        typeof(moisture_factor),
#        }(impedance_factor, viscosity_factor, moisture_factor)
#end

