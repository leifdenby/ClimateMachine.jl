# abstract type AbstractHeatModel <: BalanceLaw end

# """
#     SoilHeatModel <: BalanceLaw

# The balance law for internal energy in soil.

# """
# struct SoilHeatModel <: AbstractHeatModel end

# """
#     vars_state_conservative(soil::AbstractHeatModel, FT)

# """
# function vars_state_conservative(soil::AbstractHeatModel, FT)
#     @vars()
# end

# """
#     ConstantInternalEnergy{FT} <: AbstractHeatModel

# """
# struct ConstantInternalEnergy{FT} <: AbstractHeatModel
#     T::FT
# end

# """
#     get_temperature(m::ConstantInternalEnergy) (two methods, define it twice, one for contant energy, another for energy that comes from heat equation)

# """
# get_temperature(m::ConstantInternalEnergy) = m.T

### Soil heat model
export SoilHeatModel

"""
    SoilHeatModel{FT, SP} <: BalanceLaw

The necessary components for Heat Equation for water in soil.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SoilHeatModel{FT, SP, FiT, BCD, BCN} <: BalanceLaw # Fiθi,
    "Soil Params"
    params::SP
    "Initial conditions for temperature"
    initialT::FiT
    # "IC: volumetric ice fraction"
    # initialθ_ice::Fiθi
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
        dirichlet_bc::bc_functions = nothing,
        neumann_bc::bc_functions = nothing
    ) where {FT}

Constructor for the SoilHeatModel. Defaults imply ....
"""
function SoilHeatModel(
    ::Type{FT};
    params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
    initialT = (aux) -> FT(NaN),
    #initialθ_ice = (aux) -> FT(0.0),
    dirichlet_bc::bc_functions = nothing,
    neumann_bc::bc_functions = nothing
) where {FT}
    args = (
        params,
        initialT,
        # initialθ_ice,
        dirichlet_bc,
        neumann_bc
    )
    return SoilHeatModel{FT, typeof.(args)...}(args...)
end


"""
    vars_state_conservative(heat::SoilHeatModel, FT)

Conserved state variables (Prognostic Variables)
"""
vars_state_conservative(heat::SoilHeatModel, FT) = @vars(I::FT) #, θ_ice::FT)


"""
    vars_state_auxiliary(heat::SoilHeatModel, FT)

Names of variables required for the balance law that aren't related to derivatives of the state variables (e.g. spatial coordinates or various integrals) or those needed to solve expensive auxiliary equations (e.g., temperature via a non-linear equation solve)
"""
vars_state_auxiliary(heat::SoilHeatModel, FT) = @vars(T::FT)


"""
    vars_state_gradient(heat::SoilHeatModel, FT)

Names of the gradients of functions of the conservative state variables. Used to represent values before **and** after differentiation
"""
vars_state_gradient(heat::SoilHeatModel, FT) = @vars(I::FT)


"""
    vars_state_gradient_flux(heat::SoilHeatModel, FT)

Names of the gradient fluxes necessary to impose Neumann boundary conditions
"""
vars_state_gradient_flux(heat::SoilHeatModel, FT) = @vars(α∇I::SVector{3, FT}) ## keep in mind that rho c t is I in land doc, best to change to land document formulation


"""
    flux_first_order!(
    	land::LandModel,
    	heat::SoilHeatModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Computes and assembles non-diffusive fluxes in the model equations.
"""
function flux_first_order!(
    land::LandModel,
    heat::SoilHeatModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end


"""
    heat_init_aux!(
		land::LandModel,
		soil::SoilModel,
		heat::SoilHeatModel,
		aux::Vars,
		geom::LocalGeometry
	)

Computes and assembles non-diffusive fluxes in the model equations.
"""
function heat_init_aux!(
	land::LandModel,
	soil::SoilModel,
	heat::SoilHeatModel,
	aux::Vars,
	geom::LocalGeometry
	)
    aux.soil.heat.T = heat.initialT(aux)
end


"""
	land_nodal_update_auxiliary_state!(
    	land::LandModel,
    	soil::SoilModel,
    	heat::SoilHeatModel,
    	state::Vars,
    	aux::Vars,
    	t::Real
	)

Method of `land_nodal_update_auxiliary_state!` defining how to update the 
auxiliary variables of the soil heat balance law.
"""
function land_nodal_update_auxiliary_state!(
    land::LandModel,
    soil::SoilModel,
    heat::SoilHeatModel,
    state::Vars,
    aux::Vars,
    t::Real
)
    aux.soil.heat.T = state.soil.heat.I / heat.params.ρc
end


"""
    compute_gradient_argument!(
    	land::LandModel,
    	heat::SoilHeatModel,
        transform::Vars,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Specify how to compute the arguments to the gradients.
"""
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


"""
    compute_gradient_flux!(
        land::LandModel,
    	heat::SoilHeatModel,
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
    heat::SoilHeatModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    diffusive.soil.heat.α∇I = -heat.params.α * ∇transform.soil.heat.I
end

"""
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

Specify `F_{second_order}` for each conservative state variable
"""
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
    flux.soil.heat.I += diffusive.soil.heat.α∇I
end

