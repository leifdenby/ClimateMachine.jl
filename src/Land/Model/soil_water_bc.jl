function top_boundary_conditions!(bc::Neumann, state⁺::Vars, diff⁺::Vars,
                                  aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                                  t)
    if bc.surface_flux != nothing
        diff⁺.soil.water.κ∇h = n̂*bc.surface_flux(aux⁻, t)*aux⁻.soil.water.κ
    else
        nothing
    end
end

function top_boundary_conditions!(bc::Dirichlet, state⁺::Vars, aux⁺::Vars,
                                  state⁻::Vars, aux⁻::Vars, t
                                  )
    if bc.surface_state != nothing
        state⁺.soil.water.ϑ = bc.surface_state(aux⁻,t)
    else
        nothing
    end
end


function bottom_boundary_conditions!(bc::Neumann, state⁺::Vars, diff⁺::Vars,
                                     aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,t
                                     ) where {FT}
    if bc.bottom_flux != nothing
        diff⁺.soil.water.κ∇h = -n̂*bc.bottom_flux(aux⁻, t)*aux⁻.soil.water.κ
    else
        nothing
    end
end

function bottom_boundary_conditions!(bc::Dirichlet, state⁺::Vars, aux⁺::Vars,
                                     state⁻::Vars, aux⁻::Vars,t
                                     ) where {FT}
    if bc.bottom_state != nothing
        state⁺.soil.water.ϑ = bc.bottom_state(aux⁻, t)
    else
        nothing
    end
end



# ---------------- 4e) Boundary Conditions

# Boundary condition function
function boundary_state!(nf, land::LandModel, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...
                         )
    bc = land.soil.water.dirichlet_bc
    if bctype == 2
        top_boundary_conditions!(bc, state⁺, aux⁺, state⁻,aux⁻, t)
    elseif bctype == 1
        bottom_boundary_conditions!(bc, state⁺, aux⁺, state⁻,aux⁻, t)
    end
end

# Boundary condition function
function boundary_state!(nf, land::LandModel, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, n̂, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...
                         )
    bc = land.soil.water.neumann_bc
    if bctype == 2
        top_boundary_conditions!(bc, state⁺, diff⁺, aux⁺, n̂, state⁻, diff⁻, aux⁻,t)
    elseif bctype == 1
        bottom_boundary_conditions!(bc, state⁺, diff⁺, aux⁺, n̂, state⁻, diff⁻, aux⁻,t)
  end
end
