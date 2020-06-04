#### Turbulence model kernels

function compute_buoyancy_gradients(
    ss::SingleStack{FT, N},
    m::MixingLengthModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
) where {FT, N}
    # think how to call subdomain statistics here to get cloudy and dry values of T if you nee them
    # buoyancy gradients via chain-role
    # Alias convention:
    gm = state
    en = state
    up = state.edmf.updraft
    gm_s = source
    en_s = source
    up_s = source.edmf.updraft
    en_d = state.edmf.environment.diffusive

    _cv_d::FT     = cv_d(param_set)
    _cp_v::FT     = cp_v(param_set)
    _cp_l::FT     = cp_l(param_set)
    _cp_i::FT     = cp_i(param_set)
    _T_0::FT      = T_0(param_set)
    _e_int_i0::FT = e_int_i0(param_set)
    _grav::FT     = grav(param_set)
    _R_d::FT      = R_d(param_set)
    ϵ_v::FT       = 1 / molmass_ratio(param_set)

    cloudy, dry, cld_frac = compute_subdomain_statistics!(m, state, aux ,t, ss.statistical_model)
    ∂b∂ρ = - _grav/gm.ρ
    #                  <-------- ∂ρ∂T -------->*<----- ∂T∂e_int ---------->
    ∂ρ∂e_int_dry    = - _R_d*gm_a.p_0/(dry.R_m*dry.T*dry.T)/((1-dry.q_tot)*_cv_d+dry.q_vap *_cv_v)
    #                  <-------- ∂ρ∂T --------->*<----- ∂T∂e_int ---------->
    ∂ρ∂e_int_cloudy = - (_R_d*gm_a.p_0/(cloudy.R_m*cloudy.T*cloudy.T)/((1-cloudy.q_tot)*_cv_d+cloudy.q_vap *_cv_v+cloudy.q_liq*_cv_l+ cloudy.q_ice*_cv_i)
                       + gm_a.p_0/(cloudy.R_m*cloudy.R_m*cloudy.T)*ϵ_v*_R_d/(_cv_v*(cloudy.T-_T0)+_e_int_i0) )
    #                    <----- ∂ρ∂Rm ------->*<------- ∂Rm∂e_int ---------->

    ∂ρ∂e_int = (cld_frac * ∂ρ∂e_int_cloudy + (1-cld_frac) * ∂ρ∂e_int_dry)
    ∂ρ∂q_tot = _R_d*gm_a.p_0/(R_m*R_m*T)
    # apply chain-role
    ∂b∂z_e_int = ∂b∂ρ * ∂ρ∂e_int * ∂e_int∂z
    ∂b∂z_q_tot = ∂b∂ρ * ∂ρ∂q_tot * ∂q_tot∂z

    # add here a computation of buoyancy frequeacy based on θ_lv
    # buoyancy_freq =

    return ∂b∂z_e_int, ∂b∂z_q_tot, buoyancy_freq
end;

function gradient_Richardson_number(∂b∂z_e_int, Shear, ∂b∂z_q_tot, minval)
    Grad_Ri = min(∂b∂z_e_int/max(Shear, eps(FT)) + ∂b∂z_q_tot/max(Shear, eps(FT)) , minval)
    return Grad_Ri
end;

function turbulent_Prandtl_number(Pr_n, Grad_Ri, obukhov_length)
    if unstable(obukhov_length)
      Pr_z = Pr_n
    else
      Pr_z = Pr_n*(2*Grad_Ri/(1+(FT(53)/FT(13))*Grad_Ri -sqrt( (1+(FT(53)/FT(130))*Grad_Ri)^2 - 4*Grad_Ri)))
    end
    return Pr_z
end;

function sgs_buoyancy_freq(
    ss::SingleStack{FT, N},
    m::MixingLengthModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
    ) where {FT, N}
    # placeholder for a new buoyancy frequancy function
    return
end;

function compute_windspeed(
    ss::SingleStack{FT, N},
    m::MixingLengthModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
    ) where {FT, N}
    windspeed_min = FT(0.01) # does this needs to be exposed ?
  return max(hypot(gm.u, gm.v), windspeed_min)
end;

