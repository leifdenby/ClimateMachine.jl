include(joinpath("..","helper_funcs", "diagnose_environment.jl"))

function nondimensional_exchange_functions(
    m::AtmosModel{FT},
    entr::EntrainmentDetrainment,
    state::Vars,
    aux::Vars,
    t::Real,
    i::Int,
) where {FT, N}

    # Alias convention:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    gm_a = aux
    up_a = aux.turbconv.updraft
    en_a = aux.turbconv.environment

    # precompute vars
    _grav = FT(grav(m.param_set))
    N_upd = n_updrafts(m.turbconv)
    ρinv = 1 / gm.ρ
    up_area = up[i].ρa * ρinv
    w_up = up[i].ρau[3] / up[i].ρa
    a_en = environment_area(state, aux, N_upd)
    w_en = environment_w(state, aux, N_upd)
    en_θ_liq = environment_θ_liq(m, state, aux, N_upd)
    en_q_tot = environment_q_tot(state, aux, N_upd)

    sqrt_tke = sqrt(abs(en.ρatke) * ρinv / a_en)
    e_int = internal_energy(m, state, aux)
    ts = PhaseEquil(m.param_set, e_int, gm.ρ, gm.moisture.ρq_tot / gm.ρ)
    gm_p = air_pressure(ts)

    ts_up = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, up[i].ρaθ_liq/up[i].ρa, gm_p, up[i].ρaq_tot/up[i].ρa)
    RH_up    = relative_humidity(ts_up)

    ts_en = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, en_θ_liq, gm_p, en_q_tot)
    RH_en    = relative_humidity(ts_en)

    dw = max(w_up - w_en, 1e-4)
    db = up_a[i].buoyancy - en_a.buoyancy

    if RH_up == 1.0 || RH_en == 1.0
        c_δ = entr.c_δ
    else
        c_δ = 0.0
    end
    # compute dry and moist aux functions
    D_ε =
        entr.c_ε /
        (1 + exp(-db / dw / entr.μ_0 * (entr.χ - up_area / (up_area + a_en))))
    D_δ =
        entr.c_ε /
        (1 + exp(db / dw / entr.μ_0 * (entr.χ - up_area / (up_area + a_en))))
    M_δ = entr.c_δ * (max((RH_up^entr.β - RH_en^entr.β), 0.0))^(1 / entr.β)
    M_ε = entr.c_δ * (max((RH_en^entr.β - RH_up^entr.β), 0.0))^(1 / entr.β)
    return D_ε, D_δ, M_δ, M_ε
end;
