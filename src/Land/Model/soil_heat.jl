### Soil heat model
export SoilHeatModel

abstract type AbstractHeatModel <: BalanceLaw end

"""
    SoilHeatModel <: BalanceLaw

The balance law for internal energy in soil.

"""
struct SoilHeatModel <: AbstractHeatModel end

"""
    vars_state_conservative(soil::AbstractHeatModel, FT)

"""
function vars_state_conservative(soil::AbstractHeatModel, FT)
    @vars()
end

"""
    ConstantInternalEnergy{FT} <: AbstractHeatModel

"""
struct ConstantInternalEnergy{FT} <: AbstractHeatModel
    T::FT
end

"""
    get_temperature(m::ConstantInternalEnergy)

"""
get_temperature(m::ConstantInternalEnergy) = m.T

# ### Soil heat model
# export SoilHeatModel, SoilParamSet, Dirichlet, Neumann

# abstract type AbstractSoilParameterSet{FT <: AbstractFloat} end

# """
# Document - heat soil parameter set
# """
# Base.@kwdef struct SoilParamSet{FT} <: AbstractSoilParameterSet{FT}
#     κ_dry::FT = FT(NaN)
#     κ_sat_unfrozen::FT = FT(NaN)
#     κ_sat_frozen::FT = FT(NaN)
#     ν_om::FT = FT(NaN)
#     ν_sand::FT = FT(NaN)
#     ν_gravel::FT = FT(NaN)
#     a::FT = FT(NaN)
#     b::FT = FT(NaN)
#     porosity::FT = FT(NaN)
#     "Heat capacity"
#     ρc::FT = 1
#     "Thermal diffusivity"
#     α::FT = 0.01
# end

# abstract type bc_functions end

# Base.@kwdef struct Dirichlet{Fs, Fb} <: bc_functions
#     surface_state::Fs = nothing
#     bottom_state::Fb = nothing
# end


# Base.@kwdef struct Neumann{Fs, Fb} <: bc_functions
#     surface_flux::Fs = nothing
#     bottom_flux::Fb = nothing
# end

# """
#     SoilHeatModel{FT, SP} <: BalanceLaw

# The necessary components for Heat Equation for water in soil.

# # Fields
# $(DocStringExtensions.FIELDS)
# """
# struct SoilHeatModel{FT, SP, BCD, BCN} <: BalanceLaw
#     "Soil Params"
#     params::SP
#     "Initial conditions for temperature"
#     initialT::FT = 280
#     "Dirichlet BC structure"
#     dirichlet_bc::BCD
#     "Neumann BC structure"
#     neumann_bc::BCN
# end
# end

# #Defaults imply a constant hydraulic K = K sat
# function SoilHeatModel(
#     ::Type{FT};##we need to tell it what FT is in a non-optional way.
#     params::AbstractSoilParameterSet{FT} = SoilParamSet{FT}(),
#     initialT::FT = FT(295.15),#ideally these would be set to nothing - need to figure out how to do something like that for a float.
#     dirichlet_bc::bc_functions = nothing,
#     neumann_bc::bc_functions = nothing,
# ) where {FT}
#     args = (
#         params,
#         initialT,
#         dirichlet_bc,
#         neumann_bc,
#     )
#     return SoilHeatModel{FT, typeof.(args)...}(args...)
# end

# """
#     vars_state_conservative(heat::HeatModel, FT)

# Conserved state variables (Prognostic Variables)
# """
# vars_state_conservative(heat::HeatModel, FT) = @vars(ρcT::FT);

# """
#     vars_state_auxiliary(heat::HeatModel, FT)

# Names of variables required for the balance law that aren't related to derivatives of the state variables (e.g. spatial coordinates or various integrals) or those needed to solve expensive auxiliary equations (e.g., temperature via a non-linear equation solve)
# """
# vars_state_auxiliary(heat::HeatModel, FT) = @vars(T::FT);

# """
#     vars_state_gradient(heat::HeatModel, FT)

# Names of the gradients of functions of the conservative state variables. Used to represent values before **and** after differentiation
# """
# vars_state_gradient(heat::HeatModel, FT) = @vars(ρcT::FT);

# """
#     vars_state_gradient_flux(heat::SoilHeatModel, FT)

# Names of the gradient fluxes necessary to impose Neumann boundary conditions
# """
# vars_state_gradient_flux(heat::HeatModel, FT) = @vars(α∇ρcT::SVector{3, FT}); ## keep in mind that rho c t is I in land doc, best to change to land document formulation

# """
#     flux_first_order!(
#     	land::LandModel,
#     	soil::SoilModel,
#     	heat::SoilHeatModel,
#         flux::Grad,
#         state::Vars,
#         aux::Vars,
#         t::Real,
#     )

# Computes and assembles non-diffusive fluxes in the model equations.
# """
# function flux_first_order!(
#     land::LandModel,
#     soil::SoilModel,
#     heat::SoilHeatModel,
#     flux::Grad,
#     state::Vars,
#     aux::Vars,
#     t::Real,
#   )
# end

# """
#     heat_init_aux!(
# 		land::LandModel,
# 		soil::SoilModel,
# 		heat::SoilHeatModel,
# 		aux::Vars,
# 		geom::LocalGeometry
# 	)

# Computes and assembles non-diffusive fluxes in the model equations.
# """
# function heat_init_aux!(
# 	land::LandModel,
# 	soil::SoilModel, ## remove
# 	heat::SoilHeatModel,
# 	aux::Vars,
# 	geom::LocalGeometry
# 	)
#     aux.soil.heat.T = heat.initialT
#     #θ_ice = 0.0
#     #T = 0.0
#     # aux.soil.heat.κ = water.params.Ksat*hydraulic_conductivity(water.impedance_factor,
#     #                                                  water.viscosity_factor,
#     #                                                  water.moisture_factor,
#     #                                                  water.hydraulics,
#     #                                                  θ_ice,
#     #                                                  water.params.porosity,
#     #                                                  T,
#     #                                                  S_l)
# end

# """
# 	function land_nodal_update_auxiliary_state!(
#     	land::LandModel,
#     	soil::SoilModel,
#     	heat::SoilHeatModel,
#     	state::Vars,
#     	aux::Vars,
#     	t::Real
# 	)

# Document soil_water_nodal_update_aux! here
# """
# function land_nodal_update_auxiliary_state!(
#     land::LandModel,
#     soil::SoilModel,
#     heat::SoilHeatModel,
#     state::Vars,
#     aux::Vars,
#     t::Real
# )
#     aux.soil.heat.T = state.ρcT / m.ρc ## update this missing soil.he...
#     #θ_ice = 0.0
#     #T = 0.0
#     #S_l = effective_saturation(water.params.porosity,state.soil.water.ν)
#     #aux.soil.water.κ = water.params.Ksat*hydraulic_conductivity(water.impedance_factor,
#     #                               water.viscosity_factor,
#     #                               water.moisture_factor,
#     #                               water.hydraulics,
#     #                               θ_ice,
#     #                               water.params.porosity,
#     #                               T,
#     #                               S_l)
# end

# """
#     function compute_gradient_argument!(
#     	land::LandModel,
#     	soil::SoilModel,
#     	heat::SoilHeatModel,
#         transform::Vars,
#         state::Vars,
#         aux::Vars,
#         t::Real,
#     )

# Specify how to compute the arguments to the gradients.
# """
# function compute_gradient_argument!(
#     land::LandModel,
#     soil::SoilModel,
#     heat::SoilHeatModel,
#     transform::Vars,
#     state::Vars,
#     aux::Vars,
#     t::Real
# )
# 	transform.soil.heat.ρcT = state.soil.heat.ρcT
# end

# """
#     function compute_gradient_flux!(
#         land::LandModel,
#         soil::SoilModel,
#     	heat::SoilHeatModel,
#         diffusive::Vars,
#         ∇transform::Grad,
#         state::Vars,
#         aux::Vars,
#         t::Real,
#     )

# Specify how to compute gradient fluxes.
# """
# function compute_gradient_flux!(
#     land::LandModel,
#     soil::SoilModel,
#     heat::SoilHeatModel,
#     diffusive::Vars,
#     ∇transform::Grad,
#     state::Vars,
#     aux::Vars,
#     t::Real,
# )
#     diffusive.soil.heat.α∇ρcT = -soil.heat.α * ∇transform.soil.heat.ρcT
# end

# """
#     function flux_second_order!(
#         land::LandModel,
#         soil::SoilModel,
#     	heat::SoilHeatModel,
#         flux::Grad,
#         state::Vars,
#         diffusive::Vars,
#         hyperdiffusive::Vars,
#         aux::Vars,
#         t::Real,
#     )

# Specify `F_{second_order}` for each conservative state variable
# """
# function flux_second_order!(
#     land::LandModel,
#     soil::SoilModel,
#     heat::SoilHeatModel,
#     flux::Grad,
#     state::Vars,
#     diffusive::Vars,
#     hyperdiffusive::Vars,
#     aux::Vars,
#     t::Real,
# )
#     flux.soil.heat.ρcT += diffusive.soil.heat.α∇ρcT
# end

# # # Boundary condition function ### define in the soil instead of 
# # function boundary_state!(
# # 	nf,
# # 	land::LandModel,
# # 	state⁺::Vars,
# # 	aux⁺::Vars,
# # 	nM,
# # 	state⁻::Vars,
# # 	aux⁻::Vars,
# # 	bctype,
# # 	t,
# # 	_...)

# #   if bctype == 2
# #     # surface
# #       state⁺.soil.heat.T= land.soil.heat.surfaceT#(state⁻, aux⁻, t)
# #   elseif bctype == 1
# #     # bottom
# #     nothing
# #   end
# # end

# # # Boundary condition function
# # function boundary_state!(
# # 	nf,
# # 	land::LandModel,
# # 	state⁺::Vars,
# # 	diff⁺::Vars,
# # 	aux⁺::Vars,
# # 	n̂,
# # 	state⁻::Vars,
# # 	diff⁻::Vars,
# # 	aux⁻::Vars,
# # 	bctype,
# # 	t,
# # 	_...)
# #   if bctype == 2
# #     # surface
# #     state⁺.soil.heat.T = heat.surfaceT
# #   elseif bctype == 1
# #     # bottom
# #     diff⁺.soil.heat.α∇ρcT = n⁻ * heat.flux_bottom
# #   end
# # end

