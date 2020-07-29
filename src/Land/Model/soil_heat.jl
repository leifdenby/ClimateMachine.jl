### Soil heat model

export SoilHeatModel, PrescribedTemperatureModel

abstract type AbstractTemperatureModel <: BalanceLaw end

"""
    struct PrescribedTemperatureModel{FT} <: AbstractTemperatureModel

Model structure for a prescribed temperature model.

Document.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct PrescribedTemperatureModel{FT} <: AbstractTemperatureModel
    "Temperature"
    T::FT
end

"""
    PrescribedTemperatureModel(
        ::Type{FT};
        T = (aux,t) -> FT(0.0)
    ) where {FT}

Outer constructor for the PrescribedTemperatureModel defining default values, and
making it so changes to those defaults are supplied via keyword args.

The functions supplied by the user are point-wise evaluated and are 
evaluated in the Balance Law functions (kernels?) compute_gradient_argument,
 nodal_update, etc. whenever the prescribed temperature content variables are 
needed by the water model.
"""
function PrescribedTemperatureModel(
    ::Type{FT};
    T = (aux,t) -> FT(0.0)
) where {FT}
    args = (T)
    return PrescribedTemperatureModel{FT, typeof.(args)...}(args...)
end

"""
    SoilHeatModel{FT, SP} <: BalanceLaw

The necessary components for the Heat Equation in a soil water matrix.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilHeatModel{FT, SP, FiT, BCD, BCN} <:
        AbstractTemperatureModel
    "Soil Params"
    params::SP
    "Initial conditions for temperature"
    initialT::FiT
    "Dirichlet BC structure"
    dirichlet_bc::BCD
    "Neumann BC structure"
    neumann_bc::BCN
end

"""
    SoilHeatModel(
        ::Type{FT};
        params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
        initialT::FT = FT(NaN),
        dirichlet_bc::AbstractBoundaryFunctions = nothing,
        neumann_bc::AbstractBoundaryFunctions = nothing
    ) where {FT}

Constructor for the SoilHeatModel. Defaults imply ....
"""
function SoilHeatModel(
    ::Type{FT};
    params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
    initialT = (aux) -> FT(NaN),
    dirichlet_bc::AbstractBoundaryFunctions = nothing,
    neumann_bc::AbstractBoundaryFunctions = nothing
) where {FT}
    args = (
        params,
        initialT,
        dirichlet_bc,
        neumann_bc
    )
    return SoilHeatModel{FT, typeof.(args)...}(args...)
end


vars_state_conservative(heat::SoilHeatModel, FT) = @vars(I::FT)

vars_state_auxiliary(heat::SoilHeatModel, FT) = @vars(T::FT)
# T instead?
vars_state_gradient(heat::SoilHeatModel, FT) = @vars(I::FT)
#κ∇T instead?
vars_state_gradient_flux(heat::SoilHeatModel, FT) = @vars(α∇I::SVector{3, FT})

function flux_first_order!(
    land::LandModel,
    heat::SoilHeatModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end

function heat_init_aux!(
	land::LandModel,
	soil::SoilModel,
	heat::SoilHeatModel,
	aux::Vars,
	geom::LocalGeometry
	)
    aux.soil.heat.T = heat.initialT(aux)
end

function land_nodal_update_auxiliary_state!(
    land::LandModel,
    soil::SoilModel,
    heat::SoilHeatModel,
    state::Vars,
    aux::Vars,
    t::Real
)
    #put latent heat of fusion in - need to get initial ice content.
    aux.soil.heat.T = state.soil.heat.I / heat.params.ρc
end

function compute_gradient_argument!(
    land::LandModel,
    heat::SoilHeatModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real
)

	transform.soil.heat.I = state.soil.heat.I
end

function compute_gradient_flux!(
    land::LandModel,
    heat::SoilHeatModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    #for consistency, put the minus sign in the flux term
    diffusive.soil.heat.α∇I = heat.params.α * ∇transform.soil.heat.I
end

function flux_second_order!(
    land::LandModel,
    heat::SoilHeatModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
    # Density of liquid water (kg/m``^3``)
    _ρ_l = FT(ρ_cloud_liq(land.param_set))
    # Volum. isoboric heat capacity liquid water (J/m3/K)
    _cp_l = FT(cp_l(land.param_set) * _ρ_l)
    # Reference temperature (K)
    _T_ref = FT(T_0(land.param_set))

    # if water
    #     I_l = internal_energy_liquid_water(_cp_l, aux.soil.heat.T, _T_ref, _ρ_l)
    #     flux.soil.heat.I += - diffusive.soil.heat.α∇I - I_l * diffusive.soil.water.K∇h
    # else
        flux.soil.heat.I += - diffusive.soil.heat.α∇I
    # end

end

#Repeat for PrescribedTemperatureModel

"""
    vars_state_conservative(heat::PrescribedTemperatureModel, FT)

Conserved state variables (Prognostic Variables)
"""
vars_state_conservative(heat::PrescribedTemperatureModel, FT) = @vars()


"""
    vars_state_auxiliary(heat::PrescribedTemperatureModel, FT)

Names of variables required for the balance law that aren't related to derivatives
of the state variables (e.g. spatial coordinates or various integrals) or those
needed to solve expensive auxiliary equations (e.g., temperature via a non-linear
equation solve)
"""
vars_state_auxiliary(heat::PrescribedTemperatureModel, FT) = @vars()


"""
    vars_state_gradient(heat::PrescribedTemperatureModel, FT)

Names of the gradients of functions of the conservative state variables. Used to
represent values before **and** after differentiation
"""
vars_state_gradient(heat::PrescribedTemperatureModel, FT) = @vars()


"""
    vars_state_gradient_flux(heat::PrescribedTemperatureModel, FT)

Names of the gradient fluxes necessary to impose Neumann boundary conditions
"""
vars_state_gradient_flux(heat::PrescribedTemperatureModel, FT) = @vars()

"""
    flux_first_order!(
        land::LandModel,
        heat::PrescribedTemperatureModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Computes and assembles non-diffusive fluxes in the model equations.
"""
function flux_first_order!(
    land::LandModel,
    heat::PrescribedTemperatureModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) end


"""
    heat_init_aux!(
        land::LandModel,
        soil::SoilModel,
        heat::PrescribedTemperatureModel,
        aux::Vars,
        geom::LocalGeometry,
    )

Function defining how to initiate the auxiliary variables of the soil water balance law.
"""
function heat_init_aux!(
    land::LandModel,
    soil::SoilModel,
    heat::PrescribedTemperatureModel,
    aux::Vars,
    geom::LocalGeometry,
) end

"""
    land_nodal_update_auxiliary_state!(
        land::LandModel,
        soil::SoilModel,
        heat::PrescribedTemperatureModel
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
    heat::PrescribedTemperatureModel,
    state::Vars,
    aux::Vars,
    t::Real,
)

end

"""
    compute_gradient_argument!(
        land::LandModel,
        heat::PrescribedTemperatureModel,
        transform::Vars,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Specify how to compute the arguments to the gradients.
"""
function compute_gradient_argument!(
    land::LandModel,
    heat::PrescribedTemperatureModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)

end

"""
    compute_gradient_flux!(
        land::LandModel,
        heat::PrescribedTemperatureModel,
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
    heat::PrescribedTemperatureModel,
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
        heat::PrescribedTemperatureModel,
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
    heat::PrescribedTemperatureModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)

end

"""
    get_temperature(m::PrescribedTemperatureModel)

Returns the temperature when the heat model chosen is a user prescribed one.

This is useful for driving Richard's equation without a back reaction on temperature.
"""
get_temperature(m::PrescribedTemperatureModel) = m.T
