module Land

using CLIMAParameters
using DocStringExtensions
using LinearAlgebra, StaticArrays
using ..VariableTemplates
using ..MPIStateArrays

import ClimateMachine.DGMethods:
    BalanceLaw,
    vars_state_auxiliary,
    vars_state_conservative,
    vars_state_gradient,
    vars_state_gradient_flux,
    flux_first_order!,
    flux_second_order!,
    source!,
    boundary_state!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    init_state_auxiliary!,
    init_state_conservative!,
    update_auxiliary_state!,
    LocalGeometry,
    DGModel,
    nodal_update_auxiliary_state!


export LandModel

"""
    LandModel{PS, S, SRC} <: BalanceLaw

"""
struct LandModel{PS, S, SRC} <: BalanceLaw
    "Parameter set"
    param_set::PS
    "Soil model"
    soil::S
    "Source Terms (Problem specific source terms)"
    source::SRC
end

"""
    vars_state_conservative(land::LandModel, FT)

Conserved state variables (Prognostic Variables)
"""
function vars_state_conservative(land::LandModel, FT)
    @vars begin
        soil::vars_state_conservative(land.soil, FT)
    end
end

"""
    vars_state_auxiliary(land::LandModel, FT)
"""
function vars_state_auxiliary(land::LandModel, FT)
    @vars begin
        soil::vars_state_auxiliary(land.soil, FT)
    end
end

"""
    vars_state_gradient(land::LandModel, FT)
"""
function vars_state_gradient(land::LandModel, FT)
    @vars begin
        soil::vars_state_gradient(land.soil, FT)
    end
end

"""
    vars_state_gradient_flux(land::LandModel, FT)
"""
function vars_state_gradient_flux(land::LandModel, FT)
    @vars begin
        soil::vars_state_gradient_flux(land.soil, FT)
    end
end

function flux_first_order!(
    land::LandModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) end

function compute_gradient_argument!(
    land::LandModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)

    compute_gradient_argument!(land, land.soil, transform, state, aux, t)
end

function compute_gradient_flux!(
    land::LandModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)

    compute_gradient_flux!(
        land,
        land.soil,
        diffusive,
        ∇transform,
        state,
        aux,
        t,
    )

end

function flux_second_order!(
    land::LandModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
    flux_second_order!(
        land,
        land.soil,
        flux,
        state,
        diffusive,
        hyperdiffusive,
        aux,
        t,
    )

end

function update_auxiliary_state!(
    dg::DGModel,
    land::LandModel,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
    nodal_update_auxiliary_state!(
        land_nodal_update_auxiliary_state!,
        dg,
        m,
        Q,
        t,
        elems,
    )
end

function land_nodal_update_auxiliary_state!(
    land::LandModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    land_nodal_update_auxiliary_state!(land, land.soil, state, aux, t)
end

function source!(
    land::LandModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    land_source!(land.source, land, source, state, diffusive, aux, t, direction)
end

include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("source.jl")
include("soil_model.jl")
include("soil_heat.jl")
include("soil_water.jl")




end
