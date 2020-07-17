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
    en = state.turbulence.environment
    up = state.turbulence.updraft
    gm_a = aux
    up_a = aux.turbulence.updraft

    # precompute vars
    _grav = FT(grav(m.param_set))
    ρinv = 1 / gm.ρ
    up_area = up[i].ρa * ρinv
    w_up = up[i].ρau[3] / up[i].ρa
    en_area = 1 - sum([up[j].ρa for j in 1:N]) * ρinv
    w_en = (gm.ρu[3] - sum([up[j].ρau[3] for j in 1:N])) * ρinv
    b_up = up_a[i].buoyancy
    b_en = (
        gm_a.buoyancy -
        sum([ρinv * up[j].ρa / ρinv * up_a[j].buoyancy for j in 1:N])
    )
    sqrt_tke = sqrt(abs(en.ρatke) * ρinv / en_area)
    e_int = internal_energy(ss, state, aux)
    ts = PhaseEquil(m.param_set, e_int, gm.ρ, gm.moisture.ρq_tot / gm.ρ)
    gm_p = air_pressure(ts)

    ts_up = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, up[i].ρaθ_liq/up[i].ρa, gm_p, up[i].ρaq_tot/up[i].ρa)
    RH_up    = relative_humidity(ts_up)

    ts_en = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, en_θ_liq, gm_p, en_q_tot)
    RH_en    = relative_humidity(ts_en)

    dw = max(w_up - w_en, 1e-4)
    db = b_up - b_en

    if RH_up == 1.0 || RH_en == 1.0
        c_δ = entr.c_δ
    else
        c_δ = 0.0
    end
    # compute dry and moist aux functions
    D_ε =
        entr.c_ε /
        (1 + exp(-db / dw / entr.μ_0 * (m.χ - up_area / (up_area + en_area))))
    D_δ =
        entr.c_ε /
        (1 + exp(db / dw / entr.μ_0 * (m.χ - up_area / (up_area + en_area))))
    M_δ = entr.c_δ * (max((RH_up^m.β - RH_en^m.β), 0.0))^(1 / entr.β)
    M_ε = entr.c_δ * (max((RH_en^m.β - RH_up^m.β), 0.0))^(1 / entr.β)
    return D_ε, D_δ, M_δ, M_ε
end;
