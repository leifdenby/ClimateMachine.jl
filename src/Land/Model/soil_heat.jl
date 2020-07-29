### Soil heat model

export SoilHeatModel, PrescribedTemperatureModel

abstract type AbstractHeatModel <: AbstractSoilComponentModel end

"""
    struct PrescribedTemperatureModel{FT} <: AbstractHeatModel

Model structure for a prescribed temperature model.

Document.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct PrescribedTemperatureModel{FT} <: AbstractHeatModel
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
    SoilHeatModel{FT, SP} <: AbstractHeatModel

The necessary components for the Heat Equation in a soil water matrix.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilHeatModel{FT, SP, FiT, BCD, BCN} <:
        AbstractHeatModel
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
        params::AbstractSoilParameterFunctions{FT} = SoilParamSet{FT}(),
        initialT::FT = FT(NaN),
        dirichlet_bc::AbstractBoundaryFunctions = nothing,
        neumann_bc::AbstractBoundaryFunctions = nothing
    ) where {FT}

Constructor for the SoilHeatModel. Defaults imply ....
"""
function SoilHeatModel(
    ::Type{FT};
    params::AbstractSoilParameterFunctions{FT} = SoilParamSet{FT}(),
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

"""
    function get_temperature(
        aux::Vars,
        state::Vars,
        t::Real,
        heat::PrescribedTemperatureModel,
    )

Returns the temperature when the heat model chosen is a user prescribed one.

This is useful for driving Richard's equation without a back reaction on temperature.
"""
function get_temperature(
    aux::Vars,
    state::Vars,
    t::Real,
    heat::PrescribedTemperatureModel,
)
    T = heat.T(aux, t)
    return T
end

vars_state(heat::SoilHeatModel, st::Prognostic, FT) = @vars(I::FT)

vars_state(heat::SoilHeatModel, st::Auxiliary, FT) = @vars(T::FT)
# T instead? see comment below
vars_state(heat::SoilHeatModel, st::Gradient, FT) = @vars(I::FT)
#κ∇T instead? was modelling off of heat tutorial and also think that we should keep ∇I because I is the conserved quantity not T, it has something to do with the numerics (Brandon had suggested I use rho c t (I) for this before). i guess this isnt consistent with the how we are solving the water though
vars_state(heat::SoilHeatModel, st::GradientFlux, FT) = @vars(α∇I::SVector{3, FT})

function soil_init_aux!(
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
        flux.soil.heat.I += -diffusive.soil.heat.α∇I
    # end

end