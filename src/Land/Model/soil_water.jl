### Soil water model
export SoilWaterModel
#The model is at the end after everything else
#Everything until the model defn. should move into a separate file, like what Christian had

abstract type AbstractImpedanceFactor{FT<:AbstractFloat} end
abstract type AbstractViscosityFactor{FT<:AbstractFloat} end
abstract type AbstractMoistureFactor{FT<:AbstractFloat} end
"""
    Hydraulics model is used in the moisture factor in hydraulic conductivity and in the matric potential. The single hydraulics model choice sets both of these.
"""
abstract type AbstractHydraulicsModel{FT <: AbstractFloat} end

"""
    vanGenuchten{FT} <: AbstractHydraulicsModel{FT}
The necessary parameters for the van Genuchten hydraulic model; defaults are for Yolo light clay.
# Fields

"""
struct vanGenuchten{FT} <: AbstractHydraulicsModel{FT}
    "Exponent parameter - using in matric potential"
    n::FT
    "used in matric potential. The inverse of this carries units in the expression for matric potential (specify in inverse meters)."
    α::FT
    "Exponent parameter - determined by n, used in hydraulic conductivity"
    m::FT
    function vanGenuchten{FT}(;n::FT = FT(1.43), α::FT = FT(2.6)) where {FT}
        new(n, α, 1-1/FT(n))
    end
end

"""
    BrooksCorey{FT} <: AbstractHydraulicsModel{FT}
The necessary parameters for the Brooks and Corey hydraulic model.
Defaults are chosen to somewhat mirror the Havercamp/vG Yolo light clay hydraulic conductivity/matric potential.
# Fields

"""
Base.@kwdef struct BrooksCorey{FT} <: AbstractHydraulicsModel{FT}
    "ψ_b - used in matric potential. Units of meters."
    ψb::FT = FT(0.1656);
    "Exponent used in matric potential and hydraulic conductivity."
    m::FT = FT(0.5);
end

"""
    Haverkamp{FT} <: AbstractHydraulicsModel{FT}
The necessary parameters for the Haverkamp hydraulic model for Yolo light clay.
Note that this only is used in creating a hydraulic conductivity function, and another formulation for matric potential must be used.
# Fields

"""
struct Haverkamp{FT} <: AbstractHydraulicsModel{FT}
    "exponent in conductivity"
    k::FT
    "constant A (units of cm^k) using in conductivity. Our sim is in meters"
    A::FT
    "Exponent parameter - using in matric potential"
    n::FT
    "used in matric potential. The inverse of this carries units in the expression for matric potential (specify in inverse meters)."
    α::FT
    "Exponent parameter - determined by n, used in hydraulic conductivity"
    m::FT
    function Haverkamp{FT}(;k::FT = FT(1.77), A::FT = FT(124.6/100.0^1.77), n::FT = FT(1.43), α::FT = FT(2.6)) where {FT}
        new(k, A, n, α, 1-1/FT(n))
    end
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
    n = hm.n
    m = hm.m
    if S_l < 1
        K = sqrt(S_l)*(1-(1-S_l^(1/m))^m)^2
    else
        K = 1
    end
    return K
end


function moisture_factor(
    mm::MoistureDependent{FT},
    hm::BrooksCorey{FT},
    S_l::FT
)where {FT}
    ψb = hm.ψb
    m = hm.m

    if S_l < 1
        K = S_l^(2 * m + 3)
    else
        K = 1
    end
    return K
end


function moisture_factor(
    mm::MoistureDependent{FT},
    hm::Haverkamp{FT},
    S_l::FT,
    ψ::FT
)where {FT}
    k = hm.k
    A = hm.A
    if S_l<1
        K = A/(A+abs(ψ)^k)
    else
        K = 1
    end
    return K
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

function viscosity_factor(
    vm::TemperatureDependentViscosity{FT},
    T::FT
) where {FT}
    γ = vm.γ
    T_ref = vm.T_ref
    factor = FT(γ*(T-T_ref))
    Theta = FT(exp(factor))
    return Theta
end

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

function impedance_factor(
    imp::IceImpedance{FT},
    θ_ice::FT,
    porosity::FT
) where {FT}
    Ω = imp.Ω
    S_ice = θ_ice/porosity 
    gamma = FT(10.0^(-Ω*S_ice))
    return gamma
end

#Select models of interest. Ignore Brooks and Corey for now.
#Constant K model
function hydraulic_conductivity(
    impedance::NoImpedance{FT},
    viscosity::ConstantViscosity{FT},
    moisture::MoistureIndependent{FT};
) where {FT}
    K = FT(viscosity_factor(viscosity)*impedance_factor(impedance)*moisture_factor(moisture))
    return K
end

#Liquid model - vG
function hydraulic_conductivity(
    impedance::NoImpedance{FT},
    viscosity::ConstantViscosity{FT},
    moisture::MoistureDependent{FT},
    hydraulics::vanGenuchten{FT};
    S_l::FT
) where {FT}
    K = FT(viscosity_factor(viscosity)*
           impedance_factor(impedance)*
           moisture_factor(moisture, hydraulics, S_l))
    return K
end

#Liquid model - haverkamp
function hydraulic_conductivity(
    impedance::NoImpedance{FT},
    viscosity::ConstantViscosity{FT},
    moisture::MoistureDependent{FT},
    hydraulics::Haverkamp{FT};
    S_l::FT,
    ψ::FT
) where {FT}
    K = FT(viscosity_factor(viscosity)*
           impedance_factor(impedance)*
           moisture_factor(moisture, hydraulics, S_l, ψ))
    return K
end



#Liquid+Viscosity model - vG
function hydraulic_conductivity(
    impedance::NoImpedance{FT},
    viscosity::TemperatureDependentViscosity{FT},
    moisture::MoistureDependent{FT},
    hydraulics::vanGenuchten{FT};
    T::FT,
    S_l::FT
) where {FT}
    K = FT(viscosity_factor(viscosity, T)*
           impedance_factor(impedance)*
           moisture_factor(moisture, hydraulics, S_l))
    return K
end
   

#Liquid+ice - vG
function hydraulic_conductivity(
    impedance::IceImpedance{FT},
    viscosity::ConstantViscosity{FT},
    moisture::MoistureDependent{FT},
    hydraulics::vanGenuchten{FT};
    θ_ice::FT,
    porosity::FT,
    S_l::FT
) where {FT}
    K = FT(viscosity_factor(viscosity)*
           impedance_factor(impedance, θ_ice, porosity)*
           moisture_factor(moisture, hydraulics, S_l))
    return K
end

   

#Full model - vG
function hydraulic_conductivity(
    impedance::IceImpedance{FT},
    viscosity::TemperatureDependentViscosity{FT},
    moisture::MoistureDependent{FT},
    hydraulics::vanGenuchten{FT};
    θ_ice::FT,
    porosity::FT,
    T::FT,
    S_l::FT
) where {FT}
    K = FT(viscosity_factor(viscosity, T)*
           impedance_factor(impedance, θ_ice, porosity)*
           moisture_factor(moisture, hydraulics, S_l))
    return K
end







"""
    hydraulic_head(z,ψ)

Return the hydraulic head.

The hydraulic head is defined as the sum of vertical height and pressure head; meters.
"""
hydraulic_head(z,ψ) = z + ψ

"""
   effective_saturation(porosity::FT, θ_l::FT)

Compute the effective saturation of soil.

θ_l is defined to be zero or positive. If θ_l is negative, hydraulic functions that take it as an argument will return imaginary numbers, resulting in domain errors. However, it is possible that our current solver returns a negative θ_l due to numerical issues. Provide a warning in this case, and correct the value of θ_l so that the integration can proceed. We will remove this once the numerical issues are resolved.
"""
function effective_saturation(
        porosity::FT,
        θ_l::FT
        ) where {FT}

    if θ_l<0
        @show θ_l
        @warn("Augmented liquid fraction is negative - domain error. Artificially setting equal to zero to proceed. ")
        θ_l = 0
    end
    S_l = θ_l / porosity
    return S_l
end

"""
    pressure_head(
            model::AbstractHydraulicsModel{FT},
            porosity::FT,
            S_s::FT,
            θ_l::FT
        ) where {FT}

Determine the pressure head in both saturated and unsaturated soil.
"""
function pressure_head(
        model::AbstractHydraulicsModel{FT},
        porosity::FT,
        S_s::FT,
        θ_l::FT
        ) where {FT}
    
    S_l = effective_saturation(porosity, θ_l)
    if S_l < 1
        ψ = matric_potential(model, S_l)
    else
        ψ = (θ_l - porosity) / S_s
    end
    return ψ
end

"
    matric_potential(
            model::vanGenuchten{FT},
            S_l::FT
        ) where {FT}

Compute the van Genuchten function for matric potential.

"
function matric_potential(
        model::vanGenuchten{FT},
        S_l::FT
    ) where {FT}
    n = model.n
    m = model.m
    α = model.α

    if S_l < 1
        ψ_m = -((S_l^(-1 / m)-1) * α^(-n))^(1 / n)
    else
        ψ_m = 0
    end
    return ψ_m
end

"
    matric_potential(
            model::Haverkamp{FT},
            S_l::FT
        ) where {FT}

Compute the van Genuchten function as a proxy for the Haverkamp model matric potential (for testing purposes).

"
function matric_potential(
        model::Haverkamp{FT},
        S_l::FT
    ) where {FT}
    n = model.n
    m = mode.m
    α = model.α

    if S_l < 1
        ψ_m = -((S_l^(-1 / m)-1) * α^(-n))^(1 / n)
    else
        ψ_m = 0
    end
    return ψ_m
end

"
    matric_potential(
            model::BrooksCorey{FT},
            S_l::FT
        ) where {FT}

Compute the Brooks and Corey function for matric potential.
"
function matric_potential(
        model::BrooksCorey{FT},
        S_l::FT
    ) where {FT}
    ψb = model.ψb
    m = model.m

    if S_l <= 1
        ψ_m = -ψb*S_l^(-1/m)
    else
        ψ_m = ψb
    end
    return ψ_m
end



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
