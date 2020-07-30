function save_thermo_state_ingredients_updraft!(
    m::AtmosModel,
    moist::EquilMoist,
    state::Vars,
    aux::Vars,
)
    param_set = atmos.param_set
    FT = eltype(state)
    for i in 1:N_up
        p = aux.ref_state.p
        θ_liq = 
        q_tot = 
        ts = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, θ_liq, p, q_tot)
        aux.turbconv.updraft[i].T = air_temperature(ts)
    end
end

# function recover_thermo_state_from_ingredients_environment(
#     m::AtmosModel,
#     moist::EquilMoist,
#     state::Vars,
#     aux::Vars,
# )
#     # e_int = internal_energy(atmos, state, aux)
#     param_set = atmos.param_set
#     FT = eltype(state)

#     en_area  = environment_area(state, aux, N_up)
#     en_θ_liq = environment_θ_liq(m, state, aux, N_up)
#     en_q_tot = environment_q_tot(state, aux, N_up)

#     q_tot = environment_q_tot(state, aux, N)
#     T = aux.moisture
#     p = aux.ref_state.p

#     ρ = air_density(param_set, T, p, PhasePartition(q_tot))
#     q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
#     e_int = internal_energy(param_set, T, q)
#     return PhaseEquil{FT, typeof(param_set)}(param_set, e_int, ρ, q_tot, T)
# end


# function store_thermo_state_ingredients_environment!(
#     m::AtmosModel,
#     moist::EquilMoist,
#     state::Vars,
#     aux::Vars,
# )
#     # e_int = internal_energy(atmos, state, aux)
#     param_set = atmos.param_set
#     FT = eltype(state)

#     en_area  = environment_area(state, aux, N_up)
#     en_θ_liq = environment_θ_liq(m, state, aux, N_up)
#     en_q_tot = environment_q_tot(state, aux, N_up)
#     for i in 1:N_up
#         ts = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, en_θ_liq, gm_p, en_q_tot)
#         aux.turbconv.updraft[i].T = 
# end

# function recover_thermo_state_from_ingredients_environment(
#     m::AtmosModel,
#     moist::EquilMoist,
#     state::Vars,
#     aux::Vars,
# )
#     # e_int = internal_energy(atmos, state, aux)
#     param_set = atmos.param_set
#     FT = eltype(state)

#     en_area  = environment_area(state, aux, N_up)
#     en_θ_liq = environment_θ_liq(m, state, aux, N_up)
#     en_q_tot = environment_q_tot(state, aux, N_up)

#     q_tot = environment_q_tot(state, aux, N)
#     T = aux.moisture
#     p = aux.ref_state.p

#     ρ = air_density(param_set, T, p, PhasePartition(q_tot))
#     q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
#     e_int = internal_energy(param_set, T, q)
#     return PhaseEquil{FT, typeof(param_set)}(param_set, e_int, ρ, q_tot, T)
# end

function thermo_state_environment(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    return 1 - sum([state.turbconv.updraft[i].ρa for i in 1:N])/ state.ρ
end

function thermo_state_updraft(
    state::Vars,
    aux::Vars,
    i::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    up = state.turbconv.updraft
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft
    return a_en * en_a.buoyancy + sum([up_a[i].buoyancy*up[i].ρa*ρinv for i in 1:N])
end