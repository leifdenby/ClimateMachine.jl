"""
    top_boundary_conditions!(
        water_bc::Neumann,
        heat_bc::Neumann,
        state⁺::Vars,
        diff⁺::Vars,
        aux⁺::Vars,
        n̂,
        state⁻::Vars,
        diff⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Neumann boundary conditions for the top of the soil, if given.
"""
function top_boundary_conditions!(
    water_bc::Neumann,
    heat_bc::Neumann,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    t,
)
    if water_bc.water_surface_flux != nothing
        diff⁺.soil.water.κ∇h = n̂ * water_bc.water_surface_flux(aux⁻, t)
    else
        nothing
    end

    if heat_bc.heat_surface_flux != nothing
        diff⁺.soil.heat.α∇ρcT = n̂ * heat_bc.heat_surface_flux(aux⁻, t)
    else
        nothing
    end
end

"""
    top_boundary_conditions!(
        water_bc::Dirichlet,
        heat_bc::Dirichlet,
        state⁺::Vars,
        aux⁺::Vars,
        state⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Dirichlet boundary conditions for the top of the soil, if given.
"""
function top_boundary_conditions!(
    water_bc::Dirichlet,
    heat_bc::Dirichlet,
    state⁺::Vars,
    aux⁺::Vars,
    state⁻::Vars,
    aux⁻::Vars,
    t,
)
    if water_bc.water_surface_state != nothing
        state⁺.soil.water.ϑ = water_bc.water_surface_state(aux⁻, t)
    else
        nothing
    end

    if heat_bc.heat_surface_state != nothing
        state⁺.soil.heat.T = heat_bc.heat_surface_state(aux⁻, t)
    else
        nothing
    end
end

"""
    bottom_boundary_conditions!(
        water_bc::Neumann,
        heat_bc::Neumann,
        state⁺::Vars,
        diff⁺::Vars,
        aux⁺::Vars,
        n̂,
        state⁻::Vars,
        diff⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Neumann boundary conditions for the bottom of the soil, if given.
"""
function bottom_boundary_conditions!(
    water_bc::Neumann,
    heat_bc::Neumann,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    t,
) where {FT}

    if water_bc.water_bottom_flux != nothing
        diff⁺.soil.water.κ∇h = - n̂ * water_bc.water_bottom_flux(aux⁻, t)
    else
        nothing
    end

    if heat_bc.heat_bottom_flux != nothing
        diff⁺.soil.heat.α∇ρcT = - n̂ * heat_bc.heat_bottom_flux(aux⁻, t)
    else
        nothing
    end
end


"""
    bottom_boundary_conditions!(
        water_bc::Dirichlet,
        heat_bc::Dirichlet,
        state⁺::Vars,
        aux⁺::Vars,
        state⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Dirichlet boundary conditions for the bottom of the soil, if given.
"""
function bottom_boundary_conditions!(
    water_bc::Dirichlet,
    heat_bc::Dirichlet,
    state⁺::Vars,
    aux⁺::Vars,
    state⁻::Vars,
    aux⁻::Vars,
    t,
) where {FT}

    if water_bc.water_bottom_state != nothing
        state⁺.soil.water.ϑ = water_bc.water_bottom_state(aux⁻, t)
    else
        nothing
    end

    if heat_bc.heat_bottom_state != nothing
        state⁺.soil.heat.T = heat_bc.heat_bottom_state(aux⁻, t)
    else
        nothing
    end
end


"""
    boundary_state!(
        nf,
        land::LandModel,
        state⁺::Vars,
        aux⁺::Vars,
        nM,
        state⁻::Vars,
        aux⁻::Vars,
        bctype,
        t,
        _...,
    )

Provides boundary conditions for the balance law.
"""
function boundary_state!(
    nf,
    land::LandModel,
    state⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)

    bc_water = land.soil.water.dirichlet_bc
    bc_heat = land.soil.heat.dirichlet_bc
    if bctype == 2
        top_boundary_conditions!(bc_water, bc_heat, state⁺, aux⁺, state⁻, aux⁻, t)
    elseif bctype == 1
        bottom_boundary_conditions!(bc_water, bc_heat, state⁺, aux⁺, state⁻, aux⁻, t)
    end
end

"""
    boundary_state!(
        nf,
        land::LandModel,
        state⁺::Vars,
        diff⁺::Vars,
        aux⁺::Vars,
        n̂,
        state⁻::Vars,
        diff⁻::Vars,
        aux⁻::Vars,
        bctype,
        t,
        _...,
    )

Provides boundary conditions for the balance law.
"""
function boundary_state!(
    nf,
    land::LandModel,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
    bc_water = land.soil.water.neumann_bc
    bc_heat = land.soil.heat.neumann_bc
    if bctype == 2
        top_boundary_conditions!(
            bc_water,
            bc_heat,
            state⁺,
            diff⁺,
            aux⁺,
            n̂,
            state⁻,
            diff⁻,
            aux⁻,
            t,
        )
    elseif bctype == 1
        bottom_boundary_conditions!(
            bc_water,
            bc_heat,
            state⁺,
            diff⁺,
            aux⁺,
            n̂,
            state⁻,
            diff⁻,
            aux⁻,
            t,
        )
    end
end
