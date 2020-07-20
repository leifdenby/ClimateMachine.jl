#### Mixing length model kernels
include(joinpath("..","helper_funcs", "diagnose_environment.jl"))

function mixing_length(
    m::AtmosModel{FT, N},
    ml::MixingLengthModel,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    δ::AbstractArray{FT},
    εt::AbstractArray{FT},
) where {FT, N}

    # need to code / use the functions: obukhov_length, ustar, ϕ_m
    m.turbconv.surface
    # Alias convention:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    gm_a = aux
    en_d = diffusive.turbconv.environment
    N_upd = n_updrafts(m.turbconv)

    z = altitude(m, aux)
    _grav = FT(grav(m.param_set))
    ρinv = 1 / gm.ρ

    fill!(ml.L, 0)

    # precompute
    en_area  = environment_area(state, aux, N_upd)
    w_env    = environment_w(state, aux, N_upd)
    en_θ_liq = environment_θ_liq(m, state, aux, N_upd)
    en_q_tot = environment_q_tot(state, aux, N_upd)

    # TODO: check rank of `en_d.∇u`
    Shear = en_d.∇u[1] .^ 2 + en_d.∇u[2] .^ 2 + en_d.∇u[3] .^ 2 # consider scalar product of two vectors
    tke = en.ρatke * ρinv / en_area

    # bflux     = Nishizawa2018.compute_buoyancy_flux(m.param_set, ml.shf, ml.lhf, ml.T_b, q, ρinv)
    bflux = FT(1)
    θ_surf = m.turbconv.surface.T_surf
    # ustar = Nishizawa2018.compute_friction_velocity(m.param_set ,u_ave ,θ_suft ,flux ,Δz ,z_0 ,a ,Ψ_m_tol ,tol_abs ,iter_max)
    ustar = FT(0.28)
    # obukhov_length = Nishizawa2018.monin_obukhov_len(m.param_set, u, θ_surf, flux)
    obukhov_length = FT(-100)

    # buoyancy related functions
    ∂b∂z, Nˢ_eff = compute_buoyancy_gradients(m, state, diffusive, aux, t)
    Grad_Ri = gradient_Richardson_number(∂b∂z, Shear, FT(0.25)) # this parameter should be exposed in the model
    Pr_z = turbulent_Prandtl_number(FT(1), Grad_Ri, obukhov_length)

    # compute L1
    if Nˢ_eff > eps(FT)
        ml.L[1] = sqrt(ml.c_b * max(tke,FT(0))) / Nˢ_eff
    else
        ml.L[1] = eps(FT)
    end

    # compute L2 - law of the wall  - YAIR define tke_surf
    tke_surf = FT(1)
    if obukhov_length < eps(FT)
        ml.L[2] =
            (
                ml.κ * z / (
                    sqrt(tke_surf) / m.turbconv.surface.ustar /
                    m.turbconv.surface.ustar
                ) * ml.c_k
            ) * min((FT(1) - FT(100) * z / obukhov_length)^FT(0.2), 1 / ml.κ)
    else
        ml.L[2] =
            ml.κ * z / (
                sqrt(tke_surf) / m.turbconv.surface.ustar /
                m.turbconv.surface.ustar
            ) * ml.c_k
    end

    # compute L3 - entrainment detrainment sources
    # Production/destruction terms
    a = ml.c_m * (Shear - en_d.∇θ_liq[3] / Pr_z - en_d.∇q_tot[3] / Pr_z) * sqrt(abs(tke))
    # Dissipation term
    b = FT(0)
    # for i in 1:N
    #     a_up = up[i].ρa / gm.ρ
    #     w_up = up[i].ρau[3] / up[i].ρa
    #     b +=
    #         a_up * w_up * δ[i] / en_area *
    #         ((w_up - w_env) * (w_up - w_env) / 2 - tke) -
    #         a_up * w_up * (w_up - w_env) * εt[i] * w_env / en_area
    # end

    c_neg = ml.c_m * tke * sqrt(abs(tke))
    if abs(a) > eps(FT) && 4 * a * c_neg > -b^2
        l_entdet =
            max(-b / FT(2) / a + sqrt(b^2 + 4 * a * c_neg) / 2 / a, FT(0))
    elseif abs(a) < eps(FT) && abs(b) > eps(FT)
        l_entdet = c_neg / b
    else
        l_entdet = FT(0)
    end
    ml.L[3] = l_entdet

    frac_upper_bound = FT(0.1) # expose these in the model
    lower_bound = FT(1.5) # expose these in the model
    l_mix = lamb_smooth_minimum(ml.L, lower_bound, frac_upper_bound)
    return l_mix
end;
