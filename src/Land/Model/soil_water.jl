### Soil water model
#TODO - remove κ from aux eventually
export SoilWaterModel, SoilParamSet, Dirichlet, Neumann


abstract type AbstractSoilParameterSet{FT <: AbstractFloat} end

"""
Document - soil parameter set
"""
Base.@kwdef struct SoilParamSet{FT} <: AbstractSoilParameterSet{FT}
    porosity::FT = FT(NaN)
    Ksat::FT = FT(NaN)
    S_s::FT = FT(NaN)
end


abstract type bc_functions end

Base.@kwdef struct Dirichlet{Fs, Fb} <: bc_functions
    surface_state::Fs = nothing
    bottom_state::Fb = nothing
end
    

#scalar functions only currently
Base.@kwdef struct Neumann{Fs, Fb} <: bc_functions
    surface_flux::Fs = nothing
    bottom_flux::Fb = nothing
end


"""
    SoilWaterModel{FT, IF, VF, MF, HM, Fiϑ, BCD, BCN} <: BalanceLaw

The necessary components for Richard's Equation for water in soil.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilWaterModel{FT, SP,IF, VF, MF, HM, Fiϑ, BCD, BCN} <: BalanceLaw
    "Soil Params"
    params::SP
    "Impedance Factor - 1 or ice dependent"
    impedance_factor::IF
    "Viscosity Factor - 1 or temperature dependent"
    viscosity_factor::VF
    "Moisture Factor - 1 or moisture dependent"
    moisture_factor::MF
    "Hydraulics Model - used in matric potential and moisture factor of hydraulic conductivity."
    hydraulics::HM # vG, B/C, Haverkamp
    "IC: constant water content through profile [m3/m3]"
    initialϑ::Fiϑ
    "Dirichlet BC structure"
    dirichlet_bc::BCD
    "Neumann BC structure"
    neumann_bc::BCN
end

#Defaults imply a constant hydraulic K = K sat
function SoilWaterModel(
    ::Type{FT};##we need to tell it what FT is in a non-optional way.
    params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
    impedance_factor::AbstractImpedanceFactor{FT} = NoImpedance{FT}(),
    viscosity_factor::AbstractViscosityFactor{FT} = ConstantViscosity{FT}(),
    moisture_factor::AbstractMoistureFactor{FT} = MoistureIndependent{FT}(),
    hydraulics::AbstractHydraulicsModel{FT} = vanGenuchten{FT}(),
    initialϑ = (aux) -> FT(NaN),
    dirichlet_bc::bc_functions = nothing,
    neumann_bc::bc_functions = nothing,
) where {FT}
    args = (
        params,
        impedance_factor,
        viscosity_factor,
        moisture_factor,
        hydraulics,
        initialϑ,
        dirichlet_bc,
        neumann_bc,
    )
    return SoilWaterModel{FT, typeof.(args)...}(args...)
end


"""
    vars_state_conservative(water::SoilWaterModel, FT)

Conserved state variables (Prognostic Variables)
"""
vars_state_conservative(water::SoilWaterModel, FT) = @vars(ϑ::FT)


"""
    vars_state_auxiliary(water::SoilWaterModel, FT)

Names of variables required for the balance law that aren't related to derivatives
 of the state variables (e.g. spatial coordinates or various integrals) or those 
needed to solve expensive auxiliary equations (e.g., temperature via a non-linear 
equation solve)
"""
vars_state_auxiliary(water::SoilWaterModel, FT) = @vars(h::FT,κ::FT)


"""
    vars_state_gradient(water::SoilWaterModel, FT)

Names of the gradients of functions of the conservative state variables. Used to 
represent values before **and** after differentiation
"""
vars_state_gradient(water::SoilWaterModel, FT) = @vars(h::FT)


"""
    vars_state_gradient_flux(water::SoilWaterModel, FT)

Names of the gradient fluxes necessary to impose Neumann boundary conditions
"""
vars_state_gradient_flux(::SoilWaterModel, FT) = @vars(κ∇h::SVector{3,FT})#really, the flux is - κ∇h

"""
    flux_first_order!(
        land::LandModel,
        water::SoilWaterModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Computes and assembles non-diffusive fluxes in the model equations.
"""
function flux_first_order!(
    land::LandModel,
    water::SoilWaterModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end


"""
Document here
"""
function water_init_aux!(land::LandModel, soil::SoilModel,water::SoilWaterModel, aux::Vars, geom::LocalGeometry)
    θ_ice = 0.0
    T = 0.0
    S_l = effective_saturation(water.params.porosity,water.initialϑ(aux))
    ψ = pressure_head(water.hydraulics,
                      water.params.porosity,
                      water.params.S_s,
                      water.initialϑ(aux))
    aux.soil.water.h = hydraulic_head(aux.z,ψ)
    aux.soil.water.κ = water.params.Ksat*hydraulic_conductivity(
        water.impedance_factor,
        water.viscosity_factor,
        water.moisture_factor,
        water.hydraulics,
        θ_ice,
        water.params.porosity,
        T,
        S_l)
 end

"""
Document 
"""
function land_nodal_update_auxiliary_state!(
    land::LandModel,
    soil::SoilModel,
    water::SoilWaterModel,
    state::Vars,
    aux::Vars,
    t::Real
)
    θ_ice = 0.0
    T = 0.0
    S_l = effective_saturation(water.params.porosity,state.soil.water.ϑ)
    ψ = pressure_head(water.hydraulics,
                      water.params.porosity,
                      water.params.S_s,
                      state.soil.water.ϑ)
    aux.soil.water.h = hydraulic_head(aux.z,ψ)
    aux.soil.water.κ = water.params.Ksat*hydraulic_conductivity(
        water.impedance_factor,
        water.viscosity_factor,
        water.moisture_factor,
        water.hydraulics,
        θ_ice,
        water.params.porosity,
        T,
        S_l)
end

"""
    compute_gradient_argument!(
        land::LandModel,
        water::SoilWaterModel,
        transform::Vars,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Specify how to compute the arguments to the gradients.
"""
function compute_gradient_argument!(
    land::LandModel,
    water::SoilWaterModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real
)

    S_l = effective_saturation(water.params.porosity,state.soil.water.ϑ)
    ψ = pressure_head(water.hydraulics,
                      water.params.porosity,
                      water.params.S_s,
                      state.soil.water.ϑ)
    transform.soil.water.h = hydraulic_head(aux.z,ψ)

end

"""
    compute_gradient_flux!(
        land::LandModel,
        water::SoilWaterModel,
        diffusive::Vars,
        ∇transform::Grad,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Specify how to compute gradient fluxes.
"""
function compute_gradient_flux!(
    land::LandModel,
    water::SoilWaterModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    diffusive.soil.water.κ∇h = aux.soil.water.κ*∇transform.soil.water.h
end

"""
    flux_second_order!(
        land::LandModel,
        water::SoilWaterModel,
        flux::Grad,
        state::Vars,
        diffusive::Vars,
        hyperdiffusive::Vars,
        aux::Vars,
        t::Real,
    )

Specify the second order flux for each conservative state variable
"""
function flux_second_order!(
    land::LandModel,
    water::SoilWaterModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
    flux.soil.water.ϑ -= diffusive.soil.water.κ∇h
end

include("./soil_water_bc.jl")

