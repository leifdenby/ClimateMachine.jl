#### Entrainment-Detrainment kernels
include(joinpath("..","helper_funcs", "diagnose_environment.jl"))

function entr_detr(
    ss::AtmosModel{FT, N},
    m::EntrainmentDetrainment,
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
    en_a = aux.turbulence.environment
    up_a = aux.turbulence.updraft

    fill!(m.Λ, 0)
    # precompute vars
    _grav = FT(grav(ss.param_set))
    ρinv = 1 / gm.ρ
    up_area = up[i].ρa / gm.ρ
    a_en = environment_area(state, aux, N)
    w_up = up[i].ρau[3] / up[i].ρa
    w_en = (gm.ρu[3] - sum([up[j].ρau[3] for j in 1:N])) * ρinv/a_en
    b_up = up_a[i].buoyancy
    b_en = (
        gm_a.buoyancy -
        sum([ρinv * up[j].ρa / ρinv * up_a[j].buoyancy for j in 1:N])
        )
    en_θ_liq = environment_θ_liq(ss, state, aux, N)
    en_q_tot = environment_q_tot(state, aux, N)

    sqrt_tke = sqrt(abs(en.ρatke) * ρinv / a_en)
    ts = thermo_state(m, state, aux)
    gm_p = air_pressure(ts)
    ts_up = LiquidIcePotTempSHumEquil_given_pressure(ss.param_set, up[i].ρaθ_liq/up[i].ρa, gm_p, up[i].ρaq_tot/up[i].ρa)
    q_con_up = condensate(ts_up)
    ts_en = LiquidIcePotTempSHumEquil_given_pressure(ss.param_set, en_θ_liq, gm_p, en_q_tot)
    q_con_up = condensate(ts_up)

    dw = max(w_up - w_en, FT(1e-4))
    db = b_up - b_en

    if q_con_up * q_con_en > eps(FT)
        c_δ = m.c_δ
    else
        c_δ = 0
    end

    D_ε, D_δ, M_δ, M_ε =
        nondimensional_exchange_functions(ss, m, state, aux, t, i)

    m.Λ[1] = abs(db / dw)
    m.Λ[2] = m.c_λ * abs(db / (sqrt_tke + sqrt(eps(FT))))
    lower_bound = FT(0.1) # need to be moved ? 
    upper_bound = FT(0.0005)
    λ = lamb_smooth_minimum(m.Λ, lower_bound, upper_bound)

    # compute entrainment/detrainmnet components
    ε_trb = 2 * up_area * m.c_t * sqrt_tke / (w_up * up_area * up_a[i].updraft_top)
    ε_dyn = λ / w_up * (D_ε + M_ε)
    δ_dyn = λ / w_up * (D_δ + M_δ)
    return ε_dyn ,δ_dyn, ε_trb
end;
