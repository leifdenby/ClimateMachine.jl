#### Interface kernels

function nodal_update_auxiliary_state!(
    balance_law::BalanceLaw,
    local_state_prognostic::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
    local_state_gradient_flux::AbstractArray{FT, 1},
    t::AbstractFloat,
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_grad_flux = vars_state(balance_law, GradientFlux(), FT)
    nodal_update_auxiliary_state!(
        balance_law,
        Vars{vs_prognostic}(local_state_prognostic),
        Vars{vs_auxiliary}(local_state_auxiliary),
        Vars{vs_grad_flux}(local_state_gradient_flux),
        t,
    )
end

function flux_first_order!(
    balance_law::BalanceLaw,
    local_flux::AbstractArray{FT, 2},
    local_state_prognostic::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
    t::AbstractFloat,
    direction,
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    flux_first_order!(
        balance_law,
        Grad{vs_prognostic}(local_flux),
        Vars{vs_prognostic}(local_state_prognostic),
        Vars{vs_auxiliary}(local_state_auxiliary),
        t,
        direction,
    )
end

function flux_second_order!(
    balance_law::BalanceLaw,
    local_flux::AbstractArray{FT, 2},
    local_state_prognostic::AbstractArray{FT, 1},
    local_state_gradient_flux::AbstractArray{FT, 1},
    local_state_hyperdiffusion::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
    t::AbstractFloat,
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_grad_flux = vars_state(balance_law, GradientFlux(), FT)
    vs_hyperdiff = vars_state(balance_law, Hyperdiffusive(), FT)
    flux_second_order!(
        balance_law,
        Grad{vs_prognostic}(local_flux),
        Vars{vs_prognostic}(local_state_prognostic),
        Vars{vs_grad_flux}(local_state_gradient_flux),
        Vars{vs_hyperdiff}(local_state_hyperdiffusion),
        Vars{vs_auxiliary}(local_state_auxiliary),
        t,
    )

end

function source!(
    balance_law::BalanceLaw,
    local_source::AbstractArray{FT, 1},
    local_state_prognostic::AbstractArray{FT, 1},
    local_state_gradient_flux::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
    t::AbstractFloat,
    direction,
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_grad_flux = vars_state(balance_law, GradientFlux(), FT)

    source!(
        balance_law,
        Vars{vs_prognostic}(local_source),
        Vars{vs_prognostic}(local_state_prognostic),
        Vars{vs_grad_flux}(local_state_gradient_flux),
        Vars{vs_auxiliary}(local_state_auxiliary),
        t,
        direction,
    )

end

function compute_gradient_argument!(
    balance_law::BalanceLaw,
    local_transform::AbstractArray{FT, 1},
    local_state_prognostic::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
    t::AbstractFloat,
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_gradient = vars_state(balance_law, Gradient(), FT)
    compute_gradient_argument!(
        balance_law,
        Vars{vs_gradient}(local_transform),
        Vars{vs_prognostic}(local_state_prognostic),
        Vars{vs_auxiliary}(local_state_auxiliary),
        t,
    )

end

function compute_gradient_flux!(
    balance_law::BalanceLaw,
    local_state_gradient_flux::AbstractArray{FT, 1},
    local_transform_gradient::AbstractArray{FT, 2},
    local_state_prognostic::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
    t::AbstractFloat,
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_gradient = vars_state(balance_law, Gradient(), FT)
    vs_grad_flux = vars_state(balance_law, GradientFlux(), FT)
    compute_gradient_flux!(
        balance_law,
        Vars{vs_grad_flux}(local_state_gradient_flux),
        Grad{vs_gradient}(local_transform_gradient),
        Vars{vs_prognostic}(local_state_prognostic),
        Vars{vs_auxiliary}(local_state_auxiliary),
        t,
    )

end

function init_state_prognostic!(
    balance_law::BalanceLaw,
    l_state::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
    coords,
    args...,
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    init_state_prognostic!(
        balance_law,
        Vars{vs_prognostic}(l_state),
        Vars{vs_auxiliary}(local_state_auxiliary),
        coords,
        args...,
    )
end

function boundary_state!(
    numerical_flux,
    balance_law::BalanceLaw,
    state_prognostic⁺::AbstractArray{FT, 1},
    state_auxiliary⁺::AbstractArray{FT, 1},
    normal_vector::AbstractArray{FT, 1},
    state_prognostic⁻::AbstractArray{FT, 1},
    state_auxiliary⁻::AbstractArray{FT, 1},
    bctype::Int,
    t::AbstractFloat,
    state1⁻::AbstractArray{FT, 1},
    aux1⁻::AbstractArray{FT, 1},
) where {FT}
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    boundary_state!(
        numerical_flux,
        balance_law,
        Vars{vs_prognostic}(state_prognostic⁺),
        Vars{vs_auxiliary}(state_auxiliary⁺),
        normal_vector,
        Vars{vs_prognostic}(state_prognostic⁻),
        Vars{vs_auxiliary}(state_auxiliary⁻),
        bctype,
        t,
        Vars{vs_prognostic}(state1⁻),
        Vars{vs_auxiliary}(aux1⁻),
    )
end

function integral_load_auxiliary_state!(
    balance_law::BalanceLaw,
    local_kernel::AbstractArray{FT, 1},
    local_state_prognostic::AbstractArray{FT, 1},
    local_state_auxiliary::AbstractArray{FT, 1},
) where {FT}

    vs_up_integ = vars_state(balance_law, UpwardIntegrals(), FT)
    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    integral_load_auxiliary_state!(
        balance_law,
        Vars{vs_up_integ}(local_kernel),
        Vars{vs_prognostic}(local_state_prognostic),
        Vars{vs_auxiliary}(local_state_auxiliary),
    )

end

function integral_set_auxiliary_state!(
    balance_law::BalanceLaw,
    state_auxiliary::AbstractArray{FT, 1},
    local_kernel::AbstractArray{FT, 1},
) where {FT}
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_up_integ = vars_state(balance_law, UpwardIntegrals(), FT)
    integral_set_auxiliary_state!(
        balance_law,
        Vars{vs_auxiliary}(state_auxiliary),
        Vars{vs_up_integ}(local_kernel),
    )

end

function reverse_integral_load_auxiliary_state!(
    balance_law::BalanceLaw,
    l_V::AbstractArray{FT, 1},
    state::AbstractArray{FT, 1},
    state_auxiliary::AbstractArray{FT, 1},
) where {FT}

    vs_prognostic = vars_state(balance_law, Prognostic(), FT)
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_dn_integ = vars_state(balance_law, DownwardIntegrals(), FT)
    reverse_integral_load_auxiliary_state!(
        balance_law,
        Vars{vs_dn_integ}(l_V),
        Vars{vs_prognostic}(state),
        Vars{vs_auxiliary}(state_auxiliary),
    )

end

function reverse_integral_set_auxiliary_state!(
    balance_law::BalanceLaw,
    state_auxiliary::AbstractArray{FT, 1},
    l_V::AbstractArray{FT, 1},
) where {FT}
    vs_auxiliary = vars_state(balance_law, Auxiliary(), FT)
    vs_dn_integ = vars_state(balance_law, DownwardIntegrals(), FT)
    reverse_integral_set_auxiliary_state!(
        balance_law,
        Vars{vs_auxiliary}(state_auxiliary),
        Vars{vs_dn_integ}(l_V),
    )

end
