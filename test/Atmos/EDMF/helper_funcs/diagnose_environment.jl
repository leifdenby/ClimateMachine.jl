function environment_area(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    return 1 - sum([state.turbconv.updraft[i].ρa for i in 1:N])/ state.ρ
end

function environment_q_tot(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    return (state.moisture.ρq_tot  - sum([state.turbconv.updraft[i].ρaq_tot for i in 1:N])) / a_en*ρinv
end

function environment_θ_liq(
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    # e_int = internal_energy(m, state, aux)
    # ts = PhaseEquil(m.param_set, e_int, state.ρ, state.moisture.ρq_tot*ρinv)
    ts = thermo_state(m, state, aux)
    θ_liq = liquid_ice_pottemp(ts)
    return (θ_liq - sum([state.turbconv.updraft[i].ρaθ_liq*ρinv for i in 1:N])) / a_en
end

function environment_w(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    return (state.ρu[3] - sum([state.turbconv.updraft[i].ρaw for i in 1:N]))/a_en*ρinv
end

function grid_mean_b(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    up = state.turbconv.updraft
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft
    return a_en * en_a.buoyancy + sum([up_a[i].buoyancy*up[i].ρa*ρinv for i in 1:N])
end