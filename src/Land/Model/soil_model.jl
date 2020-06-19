#### Soil model

export SoilModel

"""
    SoilModel{W, H} <: BalanceLaw

"""
struct SoilModel{W, H} <: BalanceLaw
    "Water model"
    water::W
    "Heat model"
    heat::H
end

"""
    vars_state_conservative(soil::SoilModel, FT)

Conserved state variables (Prognostic Variables)
"""
function vars_state_conservative(soil::SoilModel, FT)
    @vars begin
        water::vars_state_conservative(soil.water, FT)
        heat::vars_state_conservative(soil.heat, FT)
    end
end

"""
    vars_state_auxiliary(soil::SoilModel, FT)

Conserved state variables (Prognostic Variables)
"""
function vars_state_auxiliary(soil::SoilModel, FT)
    @vars begin
        water::vars_state_auxiliary(soil.water, FT)
        heat::vars_state_auxiliary(soil.heat, FT)
    end
end

"""
    vars_state_gradient(soil::SoilModel, FT)

Conserved state variables (Prognostic Variables)
"""
function vars_state_gradient(soil::SoilModel, FT)
    @vars begin
        water::vars_state_gradient(soil.water, FT)
        heat::vars_state_gradient(soil.heat, FT)
    end
end

"""
    vars_state_gradient_flux(soil::SoilModel, FT)

Conserved state variables (Prognostic Variables)
"""
function vars_state_gradient_flux(soil::SoilModel, FT)
    @vars begin
        water::vars_state_gradient_flux(soil.water, FT)
        heat::vars_state_gradient_flux(soil.heat, FT)
    end
end

function flux_first_order!(
    land::LandModel,
    soil::SoilModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) end

function compute_gradient_argument!(
    land::LandModel,
    soil::SoilModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)

#   compute_gradient_argument!(
#     land,
#     soil.heat,
#     transform,
#     state,
#     aux,
#     t,
# )
#   compute_gradient_argument!(
#     land,
#     soil.water,
#     transform,
#     state,
#     aux,
#     t,
# )
end

function compute_gradient_flux!(
    land::LandModel,
    soil::SoilModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)

#   compute_gradient_flux!(
#     land,
#     soil.heat,
#     diffusive,
#     ∇transform,
#     state,
#     aux,
#     t,
# )
#   compute_gradient_flux!(
#     land,
#     soil.water,
#     diffusive,
#     ∇transform,
#     state,
#     aux,
#     t,
# )

end

function flux_second_order!(
    land::LandModel,
    soil::SoilModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
    # flux_second_order!(
    #     land,
    #     soil.heat,
    #     flux,
    #     state,
    #     diffusive,
    #     hyperdiffusive,
    #     aux,
    #     t,
    # )
    # flux_second_order!(
    #     land,
    #     soil.water,
    #     flux,
    #     state,
    #     diffusive,
    #     hyperdiffusive,
    #     aux,
    #     t,
    # )

end

function land_nodal_update_auxiliary_state!(
    land::LandModel,
    soil::SoilModel,
    state::Vars,
    aux::Vars,
    t::Real,
)

end
