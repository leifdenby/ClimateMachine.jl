"""
    land_soil_water_boundary_state!(
        nf,
        water::PrescibedWaterModel,
        diff⁺::Vars,
        state⁺::Vars,
        aux⁺::Vars,
        nM,
        state⁻::Vars,
        diff⁻::Vars,
        aux⁻::Vars,
        bctype,
        t,
        _...,
    )

Does nothing if a PrescribedWaterModel is chosen.
"""
function land_soil_water_boundary_state!(
    nf,
    water::PrescribedWaterModel,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)

end


"""
    land_soil_water_boundary_state!(
        nf,
        water::PrescibedWaterModel,
        state⁺::Vars,
        aux⁺::Vars,
        nM,
        state⁻::Vars,
        aux⁻::Vars,
        bctype,
        t,
        _...,
    )

Does nothing if a PrescribedWaterModel is chosen.
"""
function land_soil_water_boundary_state!(
    nf,
    water::PrescribedWaterModel,
    state⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)

end



"""
    land_soil_heat_boundary_state!(
        nf,
        heat::PrescibedTemperatureModel,
        state⁺::Vars,
        aux⁺::Vars,
        nM,
        state⁻::Vars,
        aux⁻::Vars,
        bctype,
        t,
        _...,
    )

Does nothing if a PrescribedTemperatureModel is chosen.
"""
function land_soil_heat_boundary_state!(
    nf,
    heat::PrescribedTemperatureModel,
    state⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)

end


"""
    land_soil_heat_boundary_state!(
        nf,
        heat::PrescribedTemperatureModel,
        diff⁺::Vars,
        state⁺::Vars,
        aux⁺::Vars,
        nM,
        state⁻::Vars,
        diff⁻::Vars,
        aux⁻::Vars,
        bctype,
        t,
        _...,
    )

Does nothing if a PrescribedTemperatureModel is chosen.
"""
function land_soil_heat_boundary_state!(
    nf,
    heat::PrescribedTemperatureModel,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)

end

"""
    land_soil_water_boundary_state!(
        nf,
        water::SoilWaterModel,
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
function land_soil_water_boundary_state!(
    nf,
    water::SoilWaterModel,
    state⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
    water_bc = water.dirichlet_bc
    if bctype == 2
        top_boundary_conditions!(water, water_bc, state⁺, aux⁺, state⁻, aux⁻, t)
    elseif bctype == 1
        bottom_boundary_conditions!(
            water,
            water_bc,
            state⁺,
            aux⁺,
            state⁻,
            aux⁻,
            t,
        )
    end
end

"""
    land_soil_water_boundary_state!(
        nf,
        water::SoilWaterModel,
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
function land_soil_water_boundary_state!(
    nf,
    water::SoilWaterModel,
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
    water_bc = water.neumann_bc
    if bctype == 2
        top_boundary_conditions!(
            water,
            water_bc,
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
            water,
            water_bc,
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

## heat
"""
    land_soil_heat_boundary_state!(
        nf,
        heat::SoilHeatModel,
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
function land_soil_heat_boundary_state!(
    nf,
    heat::SoilHeatModel,
    state⁺::Vars,
    aux⁺::Vars,
    nM,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
    heat_bc = heat.dirichlet_bc
    if bctype == 2
        top_boundary_conditions!(heat, heat_bc, state⁺, aux⁺, state⁻, aux⁻, t)
    elseif bctype == 1
        bottom_boundary_conditions!(
            heat,
            heat_bc,
            state⁺,
            aux⁺,
            state⁻,
            aux⁻,
            t,
        )
    end
end

"""
    land_soil_heat_boundary_state!(
        nf,
        heat::SoilHeatModel,
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
function land_soil_heat_boundary_state!(
    nf,
    heat::SoilHeatModel,
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
    heat_bc = heat.neumann_bc
    if bctype == 2
        top_boundary_conditions!(
            heat,
            heat_bc,
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
            heat,
            heat_bc,
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

"""
    top_boundary_conditions!(
        water::SoilWaterModel,
        bc::Neumann,
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
    water::SoilWaterModel,
    bc::Neumann,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.surface_flux != nothing
        diff⁺.soil.water.K∇h = n̂ * bc.surface_flux(aux⁻, t)
    else
        nothing
    end
end

"""
    top_boundary_conditions!(
        water::SoilWaterModel,
        bc::Dirichlet,
        state⁺::Vars,
        aux⁺::Vars,
        state⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Dirichlet boundary conditions for the top of the soil, if given.
"""
function top_boundary_conditions!(
    water::SoilWaterModel,
    bc::Dirichlet,
    state⁺::Vars,
    aux⁺::Vars,
    state⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.surface_state != nothing
        state⁺.soil.water.ϑ_l = bc.surface_state(aux⁻, t)
    else
        nothing
    end
end

"""
    bottom_boundary_conditions!(
        water::SoilWaterModel,
        bc::Neumann,
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
    water::SoilWaterModel,
    bc::Neumann,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.bottom_flux != nothing
        diff⁺.soil.water.K∇h = -n̂ * bc.bottom_flux(aux⁻, t)
    else
        nothing
    end
end


"""
    bottom_boundary_conditions!(
        water::SoilWaterModel,
        bc::Dirichlet,
        state⁺::Vars,
        aux⁺::Vars,
        state⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Dirichlet boundary conditions for the bottom of the soil, if given.
"""
function bottom_boundary_conditions!(
    water::SoilWaterModel,
    bc::Dirichlet,
    state⁺::Vars,
    aux⁺::Vars,
    state⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.bottom_state != nothing
        state⁺.soil.water.ϑ_l = bc.bottom_state(aux⁻, t)
    else
        nothing
    end
end

## Heat
"""
    top_boundary_conditions!(
        heat::SoilHeatModel,
        bc::Neumann,
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
    heat::SoilHeatModel,
    bc::Neumann,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.surface_flux != nothing
        diff⁺.soil.heat.α∇ρcT = n̂ * bc.surface_flux(aux⁻, t)
    else
        nothing
    end
end

"""
    top_boundary_conditions!(
        water::SoilWaterModel,
        bc::Dirichlet,
        state⁺::Vars,
        aux⁺::Vars,
        state⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Dirichlet boundary conditions for the top of the soil, if given.
"""
function top_boundary_conditions!(
    heat::SoilHeatModel,
    bc::Dirichlet,
    state⁺::Vars,
    aux⁺::Vars,
    state⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.surface_state != nothing
        state⁺.soil.heat.I = heat.params.ρc * bc.surface_state(aux⁻, t)
    else
        nothing
    end
end

"""
    bottom_boundary_conditions!(
        heat::SoilHeatModel,
        bc::Neumann,
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
    heat::SoilHeatModel,
    bc::Neumann,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.bottom_flux != nothing
        diff⁺.soil.heat.α∇ρcT = - n̂ * bc.bottom_flux(aux⁻, t)
    else
        nothing
    end
end


"""
    bottom_boundary_conditions!(
        heat::SoilHeatModel,
        bc::Dirichlet,
        state⁺::Vars,
        aux⁺::Vars,
        state⁻::Vars,
        aux⁻::Vars,
        t,
    )

Specify Dirichlet boundary conditions for the bottom of the soil, if given.
"""
function bottom_boundary_conditions!(
    heat::SoilHeatModel,
    bc::Dirichlet,
    state⁺::Vars,
    aux⁺::Vars,
    state⁻::Vars,
    aux⁻::Vars,
    t,
)
    if bc.bottom_state != nothing
        state⁺.soil.heat.I = heat.params.ρc * bc.bottom_state(aux⁻, t)
    else
        nothing
    end
end
