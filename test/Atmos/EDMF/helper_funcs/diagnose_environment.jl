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
    # ts = thermo_state(m, state, aux)
    # θ_liq = liquid_ice_pottemp(ts)
    θ_liq = FT(280)
    return (θ_liq - sum([state.turbconv.updraft[i].ρaθ_liq for i in 1:N])) / a_en*ρinv
end

function environment_u(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    return (state.ρu .- sum([state.turbconv.updraft[i].ρau for i in 1:N])) / a_en*ρinv
end

function environment_w(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    return (state.ρu[3] - sum([state.turbconv.updraft[i].ρau[3] for i in 1:N]))/a_en*ρinv
end