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

<<<<<<< HEAD
#the variable functions get used in our test scripts, so we need to export them.
export LandModel,
    vars_state_conservative, vars_state_auxiliary, vars_state_gradient_flux

"""
    LandModel{PS, S, SRC, IS} <: BalanceLaw
=======

export LandModel

"""
    LandModel{PS, S, SRC} <: BalanceLaw
>>>>>>> master

A BalanceLaw for land modeling.
Users may over-ride prescribed default values for each field.

# Usage

<<<<<<< HEAD
    LandModel(
        param_set,
        soil,
        source
        init_state_conservative
    )

=======
    LandModel{PS, S, SRC} <: BalanceLaw
>>>>>>> master

# Fields
$(DocStringExtensions.FIELDS)
"""
<<<<<<< HEAD
struct LandModel{PS, S, SRC, IS} <: BalanceLaw
=======
struct LandModel{PS, S, SRC} <: BalanceLaw
>>>>>>> master
    "Parameter set"
    param_set::PS
    "Soil model"
    soil::S
    "Source Terms (Problem specific source terms)"
    source::SRC
<<<<<<< HEAD
    "Initial Condition (Function to assign initial values of state variables)"
    init_state_conservative::IS
end

"""
    LandModel(
        param_set::AbstractParameterSet,
        soil::BalanceLaw;
        source::SRC = (),
        init_state_conservative::IS = nothing
    ) where {SRC, IS}

Constructor for the LandModel structure. 
"""
function LandModel(
    param_set::AbstractParameterSet,
    soil::BalanceLaw;
    source::SRC = (),
    init_state_conservative::IS = nothing,
) where {SRC, IS}
    @assert init_state_conservative ≠ nothing
    land = (param_set, soil, source, init_state_conservative)
    return LandModel{typeof.(land)...}(land...)
end



"""
=======
end

"""
>>>>>>> master
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

Names of variables required for the balance law that aren't related to 
derivatives of the state variables (e.g. spatial coordinates or various 
<<<<<<< HEAD
integrals) or those needed to solve expensive auxiliary equations (e.g., 
temperature via a non-linear equation solve)
"""
function vars_state_auxiliary(land::LandModel, FT)
    @vars begin
        z::FT
=======
integrals) or those needed to solve expensive auxiliary equations 
(e.g., temperature via a non-linear equation solve)
"""
function vars_state_auxiliary(land::LandModel, FT)
    @vars begin
>>>>>>> master
        soil::vars_state_auxiliary(land.soil, FT)
    end
end

"""
    vars_state_gradient(land::LandModel, FT)

Names of the gradients of functions of the conservative state 
<<<<<<< HEAD
variables. 

Used to represent values before **and** after differentiation
=======
variables.

Used to represent values before **and** after differentiation.
>>>>>>> master
"""
function vars_state_gradient(land::LandModel, FT)
    @vars begin
        soil::vars_state_gradient(land.soil, FT)
    end
end

"""
    vars_state_gradient_flux(land::LandModel, FT)

Names of the gradient fluxes necessary to impose Neumann boundary 
<<<<<<< HEAD
conditions
=======
conditions.
>>>>>>> master
"""
function vars_state_gradient_flux(land::LandModel, FT)
    @vars begin
        soil::vars_state_gradient_flux(land.soil, FT)
    end
end

<<<<<<< HEAD
"""
    init_state_auxiliary!(
        land::LandModel,
        aux::Vars,
        geom::LocalGeometry
    )

Initialise auxiliary variables for each LandModel subcomponent.

Store Cartesian coordinate information in `aux.z`.
"""
function init_state_auxiliary!(land::LandModel, aux::Vars, geom::LocalGeometry)
    aux.z = geom.coord[3]
    land_init_aux!(land, land.soil, aux, geom)
end
=======
>>>>>>> master

"""
    flux_first_order!(
        Land::LandModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
<<<<<<< HEAD
        t::Real,
        directions
=======
        t::Real
>>>>>>> master
    )

Computes and assembles non-diffusive fluxes in the model equations.
"""
function flux_first_order!(
    land::LandModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
<<<<<<< HEAD
    directions
=======
>>>>>>> master
) end


"""
    compute_gradient_argument!(
        land::LandModel,
        transform::Vars,
        state::Vars,
        aux::Vars,
        t::Real,
    )

Specify how to compute the arguments to the gradients.
"""
function compute_gradient_argument!(
    land::LandModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)

    compute_gradient_argument!(land, land.soil, transform, state, aux, t)
end

"""
    compute_gradient_flux!(
        land::LandModel,
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

"""
    flux_second_order!(
        land::LandModel,
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

"""
    update_auxiliary_state!(
        dg::DGModel,
        land::LandModel,
        Q::MPIStateArray,
        t::Real,
        elems::UnitRange,
    )

<<<<<<< HEAD
Perform any updates to the auxiliary variables needed at the beginning of each time-step.
=======
Perform any updates to the auxiliary variables needed at the 
beginning of each time-step.
>>>>>>> master
"""
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
<<<<<<< HEAD
        land,
=======
        m,
>>>>>>> master
        Q,
        t,
        elems,
    )
end

"""
    land_nodal_update_auxiliary_state!(
        land::LandModel,
        state::Vars,
        aux::Vars,
        t::Real,
    )

<<<<<<< HEAD
Update the auxiliary state array
=======
Update the auxiliary state array.
>>>>>>> master
"""
function land_nodal_update_auxiliary_state!(
    land::LandModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    land_nodal_update_auxiliary_state!(land, land.soil, state, aux, t)
end

"""
    source!(
        land::LandModel,
        source::Vars,
        state::Vars,
        diffusive::Vars,
        aux::Vars,
        t::Real,
        direction,n
    )
Computes (and assembles) source terms `S(Y)` in the balance law.
"""
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

<<<<<<< HEAD
"""
    init_state_conservative!(
        land::LandModel,
        state::Vars,
        aux::Vars,
        coords,
        t,
        args...)
Initialise state variables.
`args...` provides an option to include configuration data
(current use cases include problem constants, spline-interpolants)
"""
function init_state_conservative!(
    land::LandModel,
    state::Vars,
    aux::Vars,
    coords,
    t,
    args...,
)
    land.init_state_conservative(land, state, aux, coords, t, args...)
end

include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
# include("SoilHeatParameterizations.jl")
# using .SoilHeatParameterizations
include("source.jl")
include("land_bc.jl")
include("soil_model.jl")
include("soil_water.jl")
include("soil_heat.jl")
include("soil_bc.jl")
=======
include("SoilWaterParameterizations.jl")
using .SoilWaterParameterizations
include("source.jl")
include("soil_model.jl")
include("soil_heat.jl")
include("soil_water.jl")

>>>>>>> master
end # Module
