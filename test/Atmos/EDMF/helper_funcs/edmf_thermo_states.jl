# Convenience wrapper
save_subdomain_temperature!(m, state, aux) =
    save_subdomain_temperature!(m,m.moisture,state,aux)

function save_subdomain_temperature!(
    m::AtmosModel,
    moist::EquilMoist,
    state::Vars,
    aux::Vars,
)
    N_up = n_updrafts(m.turbconv)
    ts_gm = thermo_state(m, state, aux)
    p = air_pressure(ts_gm)
    θ_liq_en = liquid_ice_pottemp(ts_gm)
    q_tot_en = total_specific_humidity(ts_gm)
    for i in 1:N_up
        ρinv = state.ρ
        ρa_up = state.turbconv.updraft[i].ρa
        ρaθ_liq_up = state.turbconv.updraft[i].ρaθ_liq
        ρaq_tot_up = state.turbconv.updraft[i].ρaq_tot
        θ_liq_up = ρaθ_liq_up / ρa_up
        q_tot_up = ρaq_tot_up / ρa_up
        θ_liq_en -= ρaθ_liq_up*ρinv
        q_tot_en -= ρaq_tot_up*ρinv
        ts_up = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, θ_liq_up, p, q_tot_up)
        aux.turbconv.updraft[i].T = air_temperature(ts_up)
    end
    ts_en = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, θ_liq_en, p, q_tot_en)
    aux.turbconv.environment.T = air_temperature(ts_en)
    return nothing
end

# Convenience wrapper
thermo_state_up(m, state, aux, i_up) =
    thermo_state_up(m,m.moisture,state,aux,i_up)

function thermo_state_up(
    m::AtmosModel,
    moist::EquilMoist,
    state::Vars,
    aux::Vars,
    i_up::Int
    )
    FT = eltype(state)
    param_set = m.param_set

    ts_gm = thermo_state(m, state, aux)
    p = air_pressure(ts_gm)
    T = aux.turbconv.updraft[i_up].T
    q_tot = state.turbconv.updraft[i_up].ρaq_tot / state.turbconv.updraft[i_up].ρa
    ρ = air_density(param_set, T, p, PhasePartition(q_tot))
    q = PhasePartition_equil(param_set, T, ρ, q_tot, PhaseEquil)
    e_int = internal_energy(param_set, T, q)
    return PhaseEquil{FT, typeof(param_set)}(param_set, e_int, ρ, q.tot, T)
end

# Convenience wrapper
thermo_state_en(m, state, aux) =
    thermo_state_en(m, m.moisture, state, aux)

function thermo_state_en(
    m::AtmosModel,
    moist::EquilMoist,
    state::Vars,
    aux::Vars,
    )
    FT = eltype(state)
    param_set = m.param_set
    N_up = n_updrafts(m.turbconv)

    ts_gm = thermo_state(m, state, aux)
    p = air_pressure(ts_gm)
    T = aux.turbconv.environment.T
    ρinv = 1/state.ρ
    ρaq_tot_en = total_specific_humidity(ts_gm) - sum([state.turbconv.updraft[i].ρaq_tot for i in 1:N_up])*ρinv
    a_en = 1 - sum([state.turbconv.updraft[i].ρa for i in 1:N_up])*ρinv
    q_tot = ρaq_tot_en * ρinv / a_en
    ρ = air_density(param_set, T, p, PhasePartition(q_tot))
    q = PhasePartition_equil(param_set, T, ρ, q_tot, PhaseEquil)
    e_int = internal_energy(param_set, T, q)
    return PhaseEquil{FT, typeof(param_set)}(param_set, e_int, ρ, q.tot, T)
end
