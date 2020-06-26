### Soil water model
#TODO - remove κ from aux eventually, change ν to ϑ_l
export SoilWaterModel, SoilParamSet


abstract type AbstractSoilParameterSet{FT <: AbstractFloat} end

"""
Document - soil parameter set
"""
Base.@kwdef struct SoilParamSet{FT} <: AbstractSoilParameterSet{FT}
    porosity::FT = FT(0.0)##what is "nothing" for a object of type float?
    Ksat::FT = FT(0.0)
    S_s::FT = FT(0.0)
end



"""
    SoilWaterModel{FT, IF, VF, MF, HM, Fiν, Fsν} <: BalanceLaw

The necessary components for Richard's Equation for water in soil.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilWaterModel{FT, SP,IF, VF, MF, HM, Fiν, Fsν} <: BalanceLaw
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
    initialν::Fiν
    "Surface BC: constant water content [m3/m3]"
    surfaceν::Fsν
end

#Defaults imply a constant hydraulic K = K sat
function SoilWaterModel(
    ::Type{FT};##we need to tell it what FT is in a non-optional way.
    params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
    impedance_factor::AbstractImpedanceFactor{FT} = NoImpedance{FT}(),
    viscosity_factor::AbstractViscosityFactor{FT} = ConstantViscosity{FT}(),
    moisture_factor::AbstractMoistureFactor{FT} = MoistureIndependent{FT}(),
    hydraulics::AbstractHydraulicsModel{FT} = vanGenuchten{FT}(),
    initialν::FT = FT(0.0),#ideally these would be set to nothing - need to figure out how to do something like that for a float.
    surfaceν::FT = FT(0.0),
) where {FT}
    args = (
        params,
        impedance_factor,
        viscosity_factor,
        moisture_factor,
        hydraulics,
        initialν,
        surfaceν
    )
    return SoilWaterModel{FT, typeof.(args)...}(args...)
end


"""
    vars_state_conservative(water::SoilWaterModel, FT)

Conserved state variables (Prognostic Variables)
"""
vars_state_conservative(water::SoilWaterModel, FT) = @vars(ν::FT)


"""
    vars_state_auxiliary(water::SoilWaterModel, FT)

Names of variables required for the balance law that aren't related to derivatives of the state variables (e.g. spatial coordinates or various integrals) or those needed to solve expensive auxiliary equations (e.g., temperature via a non-linear equation solve)
"""
vars_state_auxiliary(water::SoilWaterModel, FT) = @vars(h::FT,κ::FT)


"""
    vars_state_gradient(water::SoilWaterModel, FT)

Names of the gradients of functions of the conservative state variables. Used to represent values before **and** after differentiation
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
    S_l = effective_saturation(water.params.porosity,water.initialν)
    ψ = pressure_head(water.hydraulics,
                      water.params.porosity,
                      water.params.S_s,
                      water.initialν)
    aux.soil.water.h = hydraulic_head(aux.z,ψ)
    aux.soil.water.κ = water.params.Ksat*hydraulic_conductivity(water.impedance_factor,
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
    S_l = effective_saturation(water.params.porosity,state.soil.water.ν)
    ψ = pressure_head(water.hydraulics,
                      water.params.porosity,
                      water.params.S_s,
                      state.soil.water.ν)
    aux.soil.water.h = hydraulic_head(aux.z,ψ)
    aux.soil.water.κ = water.params.Ksat*hydraulic_conductivity(water.impedance_factor,
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

    S_l = effective_saturation(water.params.porosity,state.soil.water.ν)
    ψ = pressure_head(water.hydraulics,
                      water.params.porosity,
                      water.params.S_s,
                      state.soil.water.ν)
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

Specify `F_{second_order}` for each conservative state variable
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
    flux.soil.water.ν -= diffusive.soil.water.κ∇h
end


"""
    boundary_state!(nf, land::LandModel, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)

First call of boundary_state! function
"""
function boundary_state!(nf, land::LandModel, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)
  if bctype == 2
    # surface
      state⁺.soil.water.ν= land.soil.water.surfaceν#(state⁻, aux⁻, t)
  elseif bctype == 1
    # bottom
    nothing
  end
end

"""
    boundary_state!(nf, land::LandModel, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, 
                         aux⁻::Vars,
                         bctype, t, _...)

Second call of boundary_state! function
"""
function boundary_state!(nf, land::LandModel, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
  if bctype == 2
    # surface
    state⁺.soil.water.ν = land.soil.water.surfaceν
  elseif bctype == 1
    # bottom
    diff⁺.soil.water.κ∇h = -n̂*1*aux⁻.soil.water.κ 
  end
end

