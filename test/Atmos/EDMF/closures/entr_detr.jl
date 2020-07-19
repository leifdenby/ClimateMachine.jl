#### Entrainment-Detrainment kernels
include(joinpath("..","helper_funcs", "diagnose_environment.jl"))

function entr_detr(
    m::AtmosModel{FT, N},
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
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft

    fill!(entr.Λ, 0)
    N_upd = n_updrafts(m.turbconv)
    # precompute vars
    _grav = FT(grav(m.param_set))
    ρinv = 1 / gm.ρ
    up_area = up[i].ρa / gm.ρ
    a_en = environment_area(state, aux, N_upd)
    w_en = environment_w(state, aux, N_upd)
    w_up = up[i].ρau[3] / up[i].ρa

    en_θ_liq = environment_θ_liq(m, state, aux, N_upd)
    en_q_tot = environment_q_tot(state, aux, N_upd)

    sqrt_tke = sqrt(abs(en.ρatke) * ρinv / a_en)
    # ts = thermo_state(m, state, aux)
    e_int = internal_energy(m, state, aux)
    ts = PhaseEquil(m.param_set, e_int, gm.ρ, gm.moisture.ρq_tot / gm.ρ)
    gm_p = air_pressure(ts)

    ts_up = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, up[i].ρaθ_liq/up[i].ρa, gm_p, up[i].ρaq_tot/up[i].ρa)
    q_con_up = condensate(ts_up)
    ts_en = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, en_θ_liq, gm_p, en_q_tot)
    q_con_en = condensate(ts_up)

    dw = max(w_up - w_en, FT(1e-4))
    db = up_a[i].buoyancy - en_a.buoyancy

    if q_con_up * q_con_en > eps(FT)
        c_δ = entr.c_δ
    else
        c_δ = 0
    end

    D_ε, D_δ, M_δ, M_ε =
        nondimensional_exchange_functions(m, entr, state, aux, t, i)

    entr.Λ[1] = abs(db / dw)
    entr.Λ[2] = entr.c_λ * abs(db / (sqrt_tke + sqrt(eps(FT))))
    lower_bound = FT(0.1) # need to be moved ? 
    upper_bound = FT(0.0005)
    λ = lamb_smooth_minimum(entr.Λ, lower_bound, upper_bound)

    # compute entrainment/detrainmnet components
    ε_trb = 2 * up_area * entr.c_t * sqrt_tke / (w_up * up_area * up_a[i].updraft_top)
    ε_dyn = λ / w_up * (D_ε + M_ε)
    δ_dyn = λ / w_up * (D_δ + M_δ)
    return ε_dyn ,δ_dyn, ε_trb
end;
