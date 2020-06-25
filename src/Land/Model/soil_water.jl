### Soil water model

export SoilWaterModel

"""
    SoilWaterModel{FT, IF, VF, MF, HM, Fiν, Fsν} <: BalanceLaw

The necessary components for Richard's Equation for water in soil.

# Fields
$(DocStringExtensions.FIELDS)
"""

struct SoilWaterModel{FT, IF, VF, MF, HM, Fiν, Fsν} <: BalanceLaw
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
    # "[m3/m3] constant water content in soil"
    # initialh::Fih = (aux) -> aux.z + ψ_0
end

#Defaults imply a constant hydraulic K = K sat
function SoilWaterModel(
    ::Type{FT};##we need to tell it what FT is in a non-optional way.
    impedance_factor::AbstractImpedanceFactor{FT} = NoImpedance{FT}(),
    viscosity_factor::AbstractViscosityFactor{FT} = ConstantViscosity{FT}(),
    moisture_factor::AbstractMoistureFactor{FT} = MoistureIndependent{FT}(),
    hydraulics::AbstractHydraulicsModel{FT} = vanGenuchten{FT}(),
    initialν::AbstractInitialCondition{FT} = initialCondition{FT}(),
    surfaceν::AbstractBoundaryConditions{FT} = boundaryConditions{FT}()
) where {FT}
    return SoilWaterModel{
        FT,##this is a call of the default constructor. To call the outer constructor we dont need to give the type arguments.
        typeof(impedance_factor),##here we get the concrete types
        typeof(viscosity_factor),
        typeof(moisture_factor),
        typeof(hydraulics),
        typeof(initialν),
        typeof(surfaceν)
    }(
        impedance_factor,
        viscosity_factor,
        moisture_factor,
        hydraulics,
        initialν,
        surfaceν
    )

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
vars_state_auxiliary(water::SoilWaterModel, FT) = @vars(z::FT, h::FT,κ::FT)

"""
    vars_state_gradient(water::SoilWaterModel, FT)

Names of the gradients of functions of the conservative state variables. Used to represent values before **and** after differentiation
"""
vars_state_gradient(water::SoilWaterModel, FT) = @vars(h::FT)

"""
    vars_state_gradient_flux(water::SoilWaterModel, FT)

Names of the gradient fluxes necessary to impose Neumann boundary conditions
"""
vars_state_gradient_flux(water::SoilWaterModel, FT) = @vars(∇h::SVector{3,FT})

"""
    flux_first_order!(
        water::SoilWaterModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Computes and assembles non-diffusive fluxes in the model equations.
"""
function flux_first_order!(
    water::SoilWaterModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end

"""
    function compute_gradient_argument!(
        water::SoilWaterModel,
        transform::Vars,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Specify how to compute the arguments to the gradients.
"""
function compute_gradient_argument!(
    water::SoilWaterModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    # S_l = effective_saturation(porosity,state.ν)
    # ψ = pressure_head(m.WF.matric_pot, S_l,porosity,S_s,state.ν)
    # transform.h = hydraulic_head(aux.z,ψ)
end

"""
    function compute_gradient_flux!(
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
    water::SoilWaterModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
  #diffusive.∇h = ∇transform.h
end

"""
    function flux_second_order!(
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
    water::SoilWaterModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
  )
#   flux.ν -= aux.κ * diffusive.∇h
end

"""
    function update_auxiliary_state!(
        dg::DGModel,
        water::SoilWaterModel,
        Q::MPIStateArray,
        t::Real,
        elems::UnitRange,
    )

Perform any updates to the auxiliary variables needed at the beginning of each time-step.
"""
function update_auxiliary_state!(
    dg::DGModel,
    water::SoilWaterModel,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
  nodal_update_auxiliary_state!(land_nodal_update_auxiliary_state!, dg, water, Q, t, elems)
  return true
end

"""
    function land_nodal_update_auxiliary_state!(
        water::SoilWaterModel,
        soil::SoilModel,
        land::LandModel,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Update the auxiliary state array
"""
function land_nodal_update_auxiliary_state!(
    water::SoilWaterModel,
    soil::SoilModel,
    land::LandModel,
    state::Vars,
    aux::Vars,
    t::Real)
    # S_l = effective_saturation(porosity,state.ν)
    # ψ = pressure_head(m.WF.matric_pot, S_l,porosity,S_s,state.ν)
    # aux.h = hydraulic_head(aux.z,ψ)
    # aux.κ = hydraulic_conductivity(m.WF.hydraulic_cond,K_sat,S_l,hydraulic_head(aux.z,ψ), aux.z)#, "Havercamp")
end

"""
    function water_init_aux!(
        water::SoilWaterModel,
        soil::SoilModel,
        land::LandModel,
        aux::Vars,
        geom::LocalGeometry
    )

Initialise auxiliary variables for SoilWaterModel subcomponent.
Store Cartesian coordinate information in `aux.coord`.
"""
function water_init_aux!(
    water::SoilWaterModel,
    soil::SoilModel,
    land::LandModel,
    aux::Vars,
    geom::LocalGeometry
  )

    # aux.z = geom.coord[3]
    # aux.h = m.initialh(aux)
    # aux.κ = m.initialκ(aux)
end

# """
#     function water_init_state_conservative!(
#         model::WaterModel,
#         state::Vars,
#         aux::Vars,
#         coords,
#         t::Real
#     )

# Initialise auxiliary variables for each SoilModel subcomponent.
# Store Cartesian coordinate information in `aux.coord`.
# """
# function water_init_state_conservative!(
#     model::WaterModel,
#     state::Vars,
#     aux::Vars,
#     coords,
#     t::Real
#   )
#     # state.ν = m.initialν(state, aux)
# end