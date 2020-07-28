## Soil water model
#TODO - remove κ from aux eventually
export SoilWaterModel, PrescribedWaterModel

abstract type AbstractWaterModel <: BalanceLaw end

"""
    struct PrescribedWaterModel{FT, F1, F2} <: AbstractWaterModel

Model structure for a prescribed water content model.

The user supplies functions of space and time for both `ϑ_l` and
`θ_ice`. No auxiliary or state variables are added, no PDE is solved.
The defaults are no moisture anywhere, for all time.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct PrescribedWaterModel{FT, F1, F2} <: AbstractWaterModel
    "Augmented liquid fraction"
    ϑ_l::F1
    "Volumetric fraction of ice"
    θ_ice::F2
end

"""
    PrescribedWaterModel(
        ::Type{FT};
        ϑ_l = (aux,t) -> FT(0.0),
        θ_ice = (aux,t) -> FT(0.0),
    ) where {FT}

Outer constructor for the PrescribedWaterModel defining default values, and
making it so changes to those defaults are supplied via keyword args.

The functions supplied by the user are point-wise evaluated and are
evaluated in the Balance Law functions (kernels?) compute_gradient_argument,
 nodal_update, etc. whenever the prescribed water content variables are
needed by the heat model.
"""
function PrescribedWaterModel(
    ::Type{FT};
    ϑ_l = (aux, t) -> FT(0.0),
    θ_ice = (aux, t) -> FT(0.0),
) where {FT}
    args = (ϑ_l, θ_ice)
    return PrescribedWaterModel{FT, typeof.(args)...}(args...)
end

"""
    SoilWaterModel{FT, SP, IF, VF, MF, HM, Fiϑl, Fiθi, BCD, BCN} <: AbstractWaterModel

The necessary components for Richard's Equation for water in soil.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilWaterModel{FT, SP, IF, VF, MF, HM, Fiϑl, Fiθi, BCD, BCN} <:
        AbstractWaterModel
    "Soil Params"
    params::SP
    "Impedance Factor - will be 1 or ice dependent"
    impedance_factor::IF
    "Viscosity Factor - will be 1 or temperature dependent"
    viscosity_factor::VF
    "Moisture Factor - will be 1 or moisture dependent"
    moisture_factor::MF
    "Hydraulics Model - used in matric potential and moisture factor of hydraulic conductivity."
    hydraulics::HM
    "Initial condition: augmented liquid fraction"
    initialϑ_l::Fiϑl
    "Initial condition: volumetric ice fraction"
    initialθ_ice::Fiθi
    "Dirichlet boundary condition structure"
    dirichlet_bc::BCD
    "Neumann boundary condition structure"
    neumann_bc::BCN
end

"""
    SoilWaterModel(
        ::Type{FT};
        params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
        impedance_factor::AbstractImpedanceFactor{FT} = NoImpedance{FT}(),
        viscosity_factor::AbstractViscosityFactor{FT} = ConstantViscosity{FT}(),
        moisture_factor::AbstractMoistureFactor{FT} = MoistureIndependent{FT}(),
        hydraulics::AbstractHydraulicsModel{FT} = vanGenuchten{FT}(),
        initialϑ_l = (aux) -> FT(NaN),
        initialθ_ice = (aux) -> FT(NaN),
        dirichlet_bc::AbstractBoundaryFunctions = nothing,
        neumann_bc::AbstractBoundaryFunctions = nothing,
    ) where {FT}

Constructor for the SoilWaterModel. Defaults imply a constant K = K_sat model.
"""
function SoilWaterModel(
    ::Type{FT};
    params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
    impedance_factor::AbstractImpedanceFactor{FT} = NoImpedance{FT}(),
    viscosity_factor::AbstractViscosityFactor{FT} = ConstantViscosity{FT}(),
    moisture_factor::AbstractMoistureFactor{FT} = MoistureIndependent{FT}(),
    hydraulics::AbstractHydraulicsModel{FT} = vanGenuchten{FT}(),
    initialϑ_l = (aux) -> FT(NaN),
    initialθ_ice = (aux) -> FT(0.0),
    dirichlet_bc::AbstractBoundaryFunctions = nothing,
    neumann_bc::AbstractBoundaryFunctions = nothing,
) where {FT}
    args = (
        params,
        impedance_factor,
        viscosity_factor,
        moisture_factor,
        hydraulics,
        initialϑ_l,
        initialθ_ice,
        dirichlet_bc,
        neumann_bc,
    )
    return SoilWaterModel{FT, typeof.(args)...}(args...)
end

"""
    get_water_content(
        aux::Vars,
        state::Vars,
        t::Real,
        water::SoilWaterModel,
    )

Return the moisture variables for the balance law soil water model.
"""
function get_water_content(
    aux::Vars,
    state::Vars,
    t::Real,
    water::SoilWaterModel,
)
    return state.soil.water.ϑ_l, state.soil.water.θ_ice
end

"""
    get_water_content(
        aux::Vars,
        state::Vars,
        t::Real,
        water::PrescribedWaterModel,
    )

Return the moisture variables for the prescribed soil water model.
"""
function get_water_content(
    aux::Vars,
    state::Vars,
    t::Real,
    water::PrescribedWaterModel,
)
    ϑ_l = water.ϑ_l(aux, t)
    θ_ice = water.θ_ice(aux, t)
    return ϑ_l, θ_ice
end



"""
"""
function get_diffusive_water_flux(water::SoilWaterModel, diffusive::Vars)
    return diffusive.soil.water.K∇h
end

"""
"""
function get_diffusive_water_flux(water::PrescribedWaterModel, diffusive::Vars)
    return SVector{3, FT}(0, 0, 0)
end


"""
    vars_state_conservative(water::SoilWaterModel, FT)

Conserved state variables (Prognostic Variables)
"""
vars_state_conservative(water::SoilWaterModel, FT) = @vars(ϑ_l::FT, θ_ice::FT)


"""
    vars_state_auxiliary(water::SoilWaterModel, FT)

Names of variables required for the balance law that aren't related to derivatives
 of the state variables (e.g. spatial coordinates or various integrals) or those 
needed to solve expensive auxiliary equations (e.g., temperature via a non-linear 
equation solve)
"""
vars_state_auxiliary(water::SoilWaterModel, FT) = @vars(h::FT, K::FT)


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
vars_state_gradient_flux(::SoilWaterModel, FT) = @vars(K∇h::SVector{3, FT})#really, the flux is - K∇h

"""
    flux_first_order!(
        land::LandModel,
        water::SoilWaterModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
        t::Real,
        directions
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
    directions
) end


"""
    water_init_aux!(
        land::LandModel,
        soil::SoilModel,
        water::SoilWaterModel,
        aux::Vars,
        geom::LocalGeometry,
    )

Function defining how to initiate the auxiliary variables of the soil water balance law.
"""
function water_init_aux!(
    land::LandModel,
    soil::SoilModel,
    water::SoilWaterModel,
    aux::Vars,
    geom::LocalGeometry,
)
    T = soil.heat.initialT(aux) #get_temperature(land.soil.heat)
    S_l = effective_saturation(water.params.porosity, water.initialϑ_l(aux))
    ψ = pressure_head(
        water.hydraulics,
        water.params.porosity,
        water.params.S_s,
        water.initialϑ_l(aux),
    )
    aux.soil.water.h = hydraulic_head(aux.z, ψ)
    aux.soil.water.K =
        water.params.Ksat * hydraulic_conductivity(
            water.impedance_factor,
            water.viscosity_factor,
            water.moisture_factor,
            water.hydraulics,
            water.initialθ_ice(aux),
            water.params.porosity,
            T,
            S_l,
        )
end

"""
    land_nodal_update_auxiliary_state!(
        land::LandModel,
        soil::SoilModel,
        water::SoilWaterModel,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Method of `land_nodal_update_auxiliary_state!` defining how to update the
auxiliary variables of the soil water balance law.
"""
function land_nodal_update_auxiliary_state!(
    land::LandModel,
    soil::SoilModel,
    water::SoilWaterModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    T = soil.heat.initialT(aux) #get_temperature(land.soil.heat)
    S_l = effective_saturation(water.params.porosity, state.soil.water.ϑ_l)
    ψ = pressure_head(
        water.hydraulics,
        water.params.porosity,
        water.params.S_s,
        state.soil.water.ϑ_l,
    )
    aux.soil.water.h = hydraulic_head(aux.z, ψ)
    aux.soil.water.K =
        water.params.Ksat * hydraulic_conductivity(
            water.impedance_factor,
            water.viscosity_factor,
            water.moisture_factor,
            water.hydraulics,
            state.soil.water.θ_ice,
            water.params.porosity,
            T,
            S_l,
        )
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
    t::Real,
)

    S_l = effective_saturation(water.params.porosity, state.soil.water.ϑ_l)
    ψ = pressure_head(
        water.hydraulics,
        water.params.porosity,
        water.params.S_s,
        state.soil.water.ϑ_l,
    )
    transform.soil.water.h = hydraulic_head(aux.z, ψ)

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
    diffusive.soil.water.K∇h = aux.soil.water.K * ∇transform.soil.water.h
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
    flux.soil.water.ϑ_l -= diffusive.soil.water.K∇h
end

#Repeat for PrescribedWaterModel

"""
    vars_state_conservative(water::PrescribedWaterModel, FT)

Conserved state variables (Prognostic Variables)
"""
vars_state_conservative(water::PrescribedWaterModel, FT) = @vars()


"""
    vars_state_auxiliary(water::PrescribedWaterModel, FT)

Names of variables required for the balance law that aren't related to derivatives
 of the state variables (e.g. spatial coordinates or various integrals) or those 
needed to solve expensive auxiliary equations (e.g., temperature via a non-linear 
equation solve)
"""
vars_state_auxiliary(water::PrescribedWaterModel, FT) = @vars()


"""
    vars_state_gradient(water::PrescribedWaterModel, FT)

Names of the gradients of functions of the conservative state variables. Used to 
represent values before **and** after differentiation
"""
vars_state_gradient(water::PrescribedWaterModel, FT) = @vars()


"""
    vars_state_gradient_flux(water::PrescribedWaterModel, FT)

Names of the gradient fluxes necessary to impose Neumann boundary conditions
"""
vars_state_gradient_flux(::PrescribedWaterModel, FT) = @vars()

"""
    flux_first_order!(
        land::LandModel,
        water::PrescribedWaterModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Computes and assembles non-diffusive fluxes in the model equations.
"""
function flux_first_order!(
    land::LandModel,
    water::PrescribedWaterModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) end


"""
    water_init_aux!(
        land::LandModel,
        soil::SoilModel,
        water::PrescribedWaterModel,
        aux::Vars,
        geom::LocalGeometry,
    )

Function defining how to initiate the auxiliary variables of the soil water balance law.
"""
function water_init_aux!(
    land::LandModel,
    soil::SoilModel,
    water::PrescribedWaterModel,
    aux::Vars,
    geom::LocalGeometry,
) end

"""
    land_nodal_update_auxiliary_state!(
        land::LandModel,
        soil::SoilModel,
        water::PrescribedWaterModel,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Method of `land_nodal_update_auxiliary_state!` defining how to update the 
auxiliary variables of the soil water balance law.
"""
function land_nodal_update_auxiliary_state!(
    land::LandModel,
    soil::SoilModel,
    water::PrescribedWaterModel,
    state::Vars,
    aux::Vars,
    t::Real,
)

end

"""
    compute_gradient_argument!(
        land::LandModel,
        water::PrescribedWaterModel,
        transform::Vars,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Specify how to compute the arguments to the gradients.
"""
function compute_gradient_argument!(
    land::LandModel,
    water::PrescribedWaterModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)

end

"""
    compute_gradient_flux!(
        land::LandModel,
        water::PrescribedWaterModel,
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
    water::PrescribedWaterModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)

end

"""
    flux_second_order!(
        land::LandModel,
        water::PrescribedWaterModel,
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
    water::PrescribedWaterModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)

end
