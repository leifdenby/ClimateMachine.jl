#### Pressure model kernels
include(joinpath("..","helper_funcs", "diagnose_environment.jl"))

function perturbation_pressure(
    m::AtmosModel{FT, N},
    press::PressureModel,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
    i::Int,
) where {FT, N}

    # Alias convention:
    gm = state
    en = state
    up = state.turbconv.updraft
    up_a = aux.turbconv.updraft
    up_d = diffusive.turbconv.updraft

    ρinv    = 1 / gm.ρ
    N_upd = n_updrafts(m.turbconv)
    en_area = environment_area(state, aux, N_upd)
    w_env   = environment_w(state, aux, N_upd)
    w_up    = up[i].ρau[3] / up[i].ρa

    nh_press_buoy = -up[i].ρa * up_a[i].buoyancy * press.α_b
    nh_pressure_adv = up[i].ρa * press.α_a * w_up * up_d[i].∇u[3]
    nh_pressure_drag =
        -up[i].ρa * press.α_d * (w_up - w_env) * abs(w_up - w_env) / max(up_a[i].updraft_top, FT(500)) # this parameter should be exposed in the model

    dpdz = nh_press_buoy + nh_pressure_adv + nh_pressure_drag
    dpdz_tke_i = (w_up - w_env) * dpdz

    return dpdz, dpdz_tke_i
end;
