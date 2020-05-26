# # Eddy Diffusivity- Mass Flux test

# To put this in the form of ClimateMachine's [`BalanceLaw`](@ref
# ClimateMachine.DGMethods.BalanceLaw), we'll re-write the equation as:

# "tendency"      = - div("second order flux" + "first order flux") + "non-conservative source"
# \frac{∂ F}{∂ t} = - ∇ ⋅ ( F2 + F1 )                               + S

# where F1 is the flux-componenet that has no gradient term
# where F2 is the flux-componenet that has a  gradient term

# -------------- Subdomains:
# The model has a grid mean (the dycore stats vector),i=1:N updrafts and a single environment subdomain (subscript "0")
# The grid mean is prognostic in first moment and diagnostic in second moment.
# The updrafts are prognostic in first moment and set to zero in second moment.
# The environment is diagnostic in first moment and prognostic in second moment.

# ## Equations solved:
# -------------- First Moment Equations:
#                grid mean
# ``
#     "tendency"           "second order flux"   "first order flux"                 "non-conservative source"
# \frac{∂ ρ}{∂ t}         =                         - ∇ ⋅ (ρu)
# \frac{∂ ρ u}{∂ t}       = - ∇ ⋅ (-ρaK ∇u0       ) - ∇ ⋅ (ρu u' - ρ*MF_{u} )         + S_{surface Friction}
# \frac{∂ ρ e_{int}}{∂ t} = - ∇ ⋅ (-ρaK ∇e_{int,0}) - ∇ ⋅ (u ρ_{int} - ρ*MF_{e_int} ) + S_{microphysics}
# \frac{∂ ρ q_{tot}}{∂ t} = - ∇ ⋅ (-ρaK ∇E_{tot,0}) - ∇ ⋅ (u ρ_{tot} - ρ*MF_{q_tot} ) + S_{microphysics}
# MF_ϕ = \sum{a_i * (w_i-w0)(ϕ_i-ϕ0)}_{i=1:N}
# K is the Eddy_Diffusivity, given as a function of environmental variables
# ``

#                i'th updraft equations (no second order flux)
# ``
#     "tendency"                 "first order flux"    "non-conservative sources"
# \frac{∂ ρa_i}{∂ t}           = - ∇ ⋅ (ρu_i)         + (E_{i0}           - Δ_{i0})
# \frac{∂ ρa_i u_i}{∂ t}       = - ∇ ⋅ (ρu_i u_i')    + (E_{i0}*u_0       - Δ_{i0}*u_i)       + ↑*(ρa_i*b - a_i\frac{∂p^†}{∂z})
# \frac{∂ ρa_i e_{int,i}}{∂ t} = - ∇ ⋅ (ρu*e_{int,i}) + (E_{i0}*e_{int,0} - Δ_{i0}*e_{int,i}) + ρS_{int,i}
# \frac{∂ ρa_i q_{tot,i}}{∂ t} = - ∇ ⋅ (ρu*q_{tot,i}) + (E_{i0}*q_{tot,0} - Δ_{i0}*q_{tot,i}) + ρS_{tot,i}
# b = 0.01*(e_{int,i} - e_{int})/e_{int}
#
#                environment equations first moment
#
# a0 = 1-sum{a_i}{i=1:N}
# u0 = (⟨u⟩-sum{a_i*u_i}{i=1:N})/a0
# E_int0 = (⟨E_int⟩-sum{a_i*E_int_i}{i=1:N})/a0
# q_tot0 = (⟨q_tot⟩-sum{a_i*q_tot_i}{i=1:N})/a0
#
#                environment equations second moment
# ``
#     "tendency"           "second order flux"       "first order flux"  "non-conservative source"
# \frac{∂ ρa_0ϕ'ψ'}{∂ t} =  - ∇ ⋅ (-ρa_0⋅K⋅∇ϕ'ψ')  - ∇ ⋅ (u ρa_0⋅ϕ'ψ')   + 2ρa_0⋅K(∂_z⋅ϕ)(∂_z⋅ψ)  + (E_{i0}*ϕ'ψ' - Δ_{i0}*ϕ'ψ') + ρa_0⋅D_{ϕ'ψ',0} + ρa_0⋅S_{ϕ'ψ',0}
# ``

# --------------------- Ideal gas law and subdomain density
# ``
# T_i, q_l  = saturation adjustment(e_int, q_tot)
# TempShamEquil(e_int,q_tot,p)
# ρ_i = <p>/R_{m,i} * T_i
# b = -g(ρ_i-ρ_h)<ρ>

# where
#  - `t`        is time
#  - `z`        is height
#  - `ρ`        is the density
#  - `u`        is the 3D velocity vector
#  - `e_int`    is the internal energy
#  - `q_tot`    is the total specific humidity
#  - `K`        is the eddy diffusivity
#  - `↑`        is the upwards pointing unit vector
#  - `b`        is the buoyancy
#  - `E_{i0}`   is the entrainment rate from the enviroment into i
#  - `Δ_{i0}`   is the detrainment rate from i to the enviroment
#  - `ϕ'ψ'`     is a shorthand for \overline{ϕ'ψ'}_0 the enviromental covariance of ϕ and ψ
#  - `D`        is a covariance dissipation
#  - `S_{ϕ,i}`  is a source of ϕ in the i'th subdomain
#  - `∂_z`      is the vertical partial derivative

# --------------------- Initial Conditions
# Initial conditions are given for all variables in the grid mean, and subdomain variables assume their grid mean values
# ``
# ------- grid mean:
# ρ = hydrostatic reference state - need to compute that state
# ρu = 0
# ρe_int = convert from input profiles
# ρq_tot = convert from input profiles
# ------- updrafts:
# ρa_i = 0.1/N
# ρau = 0
# ρae_int = a*gm.ρe_int
# ρaq_tot = a*gm.ρq_tot
# ------- environment:
# cld_frac = 0.0
# `ϕ'ψ'` = initial covariance profile
# TKE = initial TKE profile
# ``

# --------------------- Boundary Conditions
#           grid mean
# ``
# surface: ρ =
# z_min: ρu = 0
# z_min: ρe_int = 300*cp ; cp=1000
# z_min: ρq_tot = 0.0
# ``

#           i'th updraft
# ``
# z_min: ρ = 1
# z_min: ρu = 0
# z_min: ρe_int = 302*cp ; cp=1000
# z_min: ρq_tot = 0.0
# ``


#### EDMF model kernels

import ClimateMachine.TurbulenceConvection:
    vars_state_auxiliary,
    vars_state_conservative,
    vars_state_gradient,
    vars_state_gradient_flux,
    init_state_auxiliary!,
    update_auxiliary_state!,
    flux_first_order!,
    flux_second_order!,
    turbconv_sources,
    turbconv_boundary_state!,
    turbconv_normal_boundary_flux_second_order!,
    compute_gradient_argument!,
    compute_gradient_flux!


include(joinpath("helper_funcs", "nondimensional_exchange_functions.jl"))
include(joinpath("helper_funcs", "lamb_smooth_minimum.jl"))
include(joinpath("helper_funcs", "subdomain_statistics.jl"))
include(joinpath("closures", "entr_detr.jl"))
include(joinpath("closures", "pressure.jl"))
include(joinpath("closures", "mixing_length.jl"))
include(joinpath("closures", "turbulence_functions.jl"))
include(joinpath("closures", "surface_functions.jl"))
# include(joinpath("closures", "micro_phys.jl"))

function vars_state_auxiliary(m::NTuple{N, Updraft}, FT) where {N}
    return Tuple{ntuple(i -> vars_state_auxiliary(m[i], FT), N)...}
end

function vars_state_auxiliary(::Updraft, FT)
    @vars(
        buoyancy::FT,
        updraft_top::FT,# YAIR - DO I NEED THIS
    )
end

function vars_state_auxiliary(::Environment, FT)
    @vars(
        cld_frac::FT,# YAIR - DO I NEED THIS
    )
end

function vars_state_auxiliary(m::EDMF, FT)
    @vars(
        environment::vars_state_auxiliary(m.environment, FT),
        updraft::vars_state_auxiliary(m.updraft, FT)
    )
end

function vars_state_conservative(::Updraft, FT)
    @vars(ρa::FT, ρau::SVector{3, FT}, ρae::FT, ρaq_tot::FT,)
end

function vars_state_conservative(::Environment, FT)
    @vars(ρatke::FT, ρae_int_cv::FT, ρaq_tot_cv::FT, ρae_int_q_tot_cv::FT,)
end

function vars_state_conservative(m::NTuple{N, Updraft}, FT) where {N}
    return Tuple{ntuple(i -> vars_state_conservative(m[i], FT), N)...}
end


function vars_state_conservative(m::EDMF, FT)
    @vars(
        environment::vars_state_conservative(m.environment, FT),
        updraft::vars_state_conservative(m.updraft, FT)
    )
end

function vars_state_gradient(::Updraft, FT)
    @vars(u::SVector{3, FT},)
end

function vars_state_gradient(::Environment, FT)
    @vars(
        e::FT,
        e_int::FT,
        q_tot::FT,
        u::SVector{3, FT},
        tke::FT,
        e_int_cv::FT,
        q_tot_cv::FT,
        e_int_q_tot_cv::FT,
        θ_ρ::FT,
    )
end


function vars_state_gradient(m::NTuple{N, Updraft}, FT) where {N}
    return Tuple{ntuple(i -> vars_state_gradient(m[i], FT), N)...}
end


function vars_state_gradient(m::EDMF, FT)
    @vars(
        environment::vars_state_gradient(m.environment, FT),
        updraft::vars_state_gradient(m.updraft, FT)
    )
end


function vars_state_gradient_flux(m::NTuple{N, Updraft}, FT) where {N}
    return Tuple{ntuple(i -> vars_state_gradient_flux(m[i], FT), N)...}
end

function vars_state_gradient_flux(::Updraft, FT)
    @vars(∇u::SMatrix{3, 3, FT, 9},)
end

function vars_state_gradient_flux(::Environment, FT)
    @vars(
        ∇e::SVector{3, FT},
        ∇e_int::SVector{3, FT},
        ∇q_tot::SVector{3, FT},
        ∇u::SMatrix{3, 3, FT, 9},
        ∇tke::SVector{3, FT},
        ∇e_int_cv::SVector{3, FT},
        ∇q_tot_cv::SVector{3, FT},
        ∇e_int_q_tot_cv::SVector{3, FT},
        ∇θ_ρ::SVector{3, FT},# used in a diagnostic equation for the mixing length
    )

end

function vars_state_gradient_flux(m::EDMF, FT)
    @vars(
        environment::vars_state_gradient_flux(m.environment, FT),
        updraft::vars_state_gradient_flux(m.updraft, FT)
    )
end

# Specify the initial values in `aux::Vars`, which are available in
# `init_state_conservative!`. Note that
# - this method is only called at `t=0`
# - `aux.z` and `aux.T` are available here because we've specified `z` and `T`
# in `vars_state_auxiliary`
function atmos_init_aux!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    aux::Vars,
    geom::LocalGeometry,
) where {FT, N}
    # println("-------- atmos_init_aux!")
    # Aliases:
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft

    # en_a.buoyancy = eps(FT)
    en_a.cld_frac = eps(FT)

    for i in 1:N
        up_a[i].buoyancy = eps(FT)
        up_a[i].updraft_top = eps(FT)
    end
    en_a.cld_frac = eps(FT)

end;

# Specify the initial values in `state::Vars`. Note that
# - this method is only called at `t=0`
function init_state_conservative!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
) where {FT, N}
    # println("-------- atmos_init_aux!")

    # Aliases:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft

    # gm.ρ = aux.ref_state.ρ # added at end

    # GCM setting - Initialize the grid mean profiles of prognostic variables (ρ,e_int,q_tot,u,v,w)
    z = altitude(m, aux)
    gm.ρ = FT(1)
    gm.ρu = SVector(0, 0, 0)
    gm.ρe = gm.ρ * 300000 # TODO: Revisit this initialization (used to be gm.ρe_int)
    gm.moisture.ρq_tot = gm.ρ * eps(FT)

    # SCM setting - need to have separate cases coded and called from a folder - see what LES does
    # a moist_thermo state is used here to convert the input θ,q_tot to e_int, q_tot profile

    a_up = FT(0.1)
    for i in 1:N
        up[i].ρa = gm.ρ * a_up
        up[i].ρau = gm.ρu * a_up
        up[i].ρae = gm.ρe * a_up
        up[i].ρaq_tot = gm.moisture.ρq_tot * a_up
    end

    # initialize environment covariance
    en.ρae_int_cv = eps(FT)
    en.ρaq_tot_cv = eps(FT)
    en.ρae_int_q_tot_cv = eps(FT)

end;

# The remaining methods, defined in this section, are called at every
# time-step in the solver by the [`BalanceLaw`](@ref
# ClimateMachine.DGMethods.BalanceLaw) framework.

include("compute_updraft_top.jl")

# Overload `update_auxiliary_state!` to call `single_stack_nodal_update_aux!`, or
# any other auxiliary methods
function update_auxiliary_state!(
    dg::DGModel,
    turbconv::EDMF,
    m::AtmosModel{FT},
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
) where {FT}
    # println("-------- update_auxiliary_state!")
    turbconv = m.turbconv

    # ----------
    # turbconv.updraft[i].updraft_top should be
    # defined after this method call, but it's
    # not complete yet...
    compute_updraft_top!(dg, m, turbconv, Q, t, elems)
    # ----------

    nodal_update_auxiliary_state!(
        edmf_stack_nodal_update_aux!,
        dg,
        m,
        Q,
        t,
        elems,
    )
end;

# Compute/update all auxiliary variables at each node. Note that
# - `aux.T` is available here because we've specified `T` in
# `vars_state_auxiliary`
function edmf_stack_nodal_update_aux!(
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}
    # println("-------- edmf_stack_nodal_update_aux!")
    turbconv = m.turbconv
    N = n_updrafts(m.turbconv)

    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft

    #   -------------  Compute buoyancies of subdomains
    ρinv = 1 / gm.ρ
    _grav::FT = grav(m.param_set)
    # b_upds = 0
    # a_upds = 0
    for i in 1:N
        # override ----------------
        # e_int_up = internal_energy(gm.ρ, up[i].ρae/up[i].ρa, up[i].ρau/up[i].ρa, _grav * gm_a.z)
        # ts = PhaseEquil(m.param_set ,e_int_up, gm.ρ, up[i].ρaq_tot/up[i].ρa)
        # ρ_i = #air_density(ts)
        # override ----------------
        ρ_i = aux.ref_state.ρ
        up_a[i].buoyancy = -_grav * (ρ_i - aux.ref_state.ρ) * ρinv
    end
    # compute the buoyancy of the environment
    en_area = 1 - sum([up[i].ρa for i in 1:N]) * ρinv
    env_e = (gm.ρe - sum([up[i].ρae * ρinv for i in 1:N])) / en_area
    env_q_tot =
        (gm.moisture.ρq_tot - sum([up[i].ρaq_tot * ρinv for i in 1:N])) / en_area

    # override ----------------
    # ts = PhaseEquil(m.param_set ,env_e, gm.ρ, env_q_tot)
    # env_ρ = air_density(ts)
    # env_q_liq = PhasePartition(ts).liq
    # env_q_ice = PhasePartition(ts).ice

    env_ρ = aux.ref_state.ρ
    env_q_liq = FT(0)
    env_q_ice = FT(0)
    # override ----------------

    b_env = -_grav * (env_ρ - aux.ref_state.ρ) * ρinv
    b_gm = en_area * b_env + sum([up_a[i].buoyancy for i in 1:N])
    # subtract the grid mean
    for i in 1:N
        up_a[i].buoyancy -= b_gm
    end
end;

# Since we have second-order fluxes, we must tell `ClimateMachine` to compute
# the gradient of `ρcT`. Here, we specify how `ρcT` is computed. Note that
#  - `transform.ρcT` is available here because we've specified `ρcT` in
#  `vars_state_gradient`
function compute_gradient_argument!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT, N}
    # println("-------- compute_gradient_argument!")
    # Aliases:
    up_t = transform.turbconv.updraft
    en_t = transform.turbconv.environment
    gm = state
    up = state.turbconv.updraft
    en = state.turbconv.environment

    for i in 1:N
        up_t[i].u = up[i].ρau / up[i].ρa
    end
    _grav::FT = grav(m.param_set)
    ts = thermo_state(m, state, aux)

    ρinv = 1 / gm.ρ
    en_area = 1 - sum([up[i].ρa for i in 1:N]) * ρinv

    en_ρe = (gm.ρe - sum([up[j].ρae for j in 1:N])) / en_area
    en_ρu = (gm.ρu .- sum([up[j].ρae for j in 1:N])) / en_area
    e_pot = gravitational_potential(m, aux)
    en_e_int = internal_energy(gm.ρ, en_ρe, en_ρu, e_pot)

    # populate gradient arguments
    en_t.e = en_ρe * ρinv
    en_t.e_int = en_e_int
    en_t.q_tot =
        (gm.moisture.ρq_tot - sum([up[i].ρaq_tot for i in 1:N])) / (en_area * gm.ρ)
    en_t.u = (gm.ρu .- sum([up[i].ρau for i in 1:N])) / (en_area * gm.ρ)

    en_t.tke = en.ρatke / (en_area * gm.ρ)
    en_t.e_int_cv = en.ρae_int_cv / (en_area * gm.ρ)
    en_t.q_tot_cv = en.ρaq_tot_cv / (en_area * gm.ρ)
    en_t.e_int_q_tot_cv = en.ρae_int_q_tot_cv / (en_area * gm.ρ)

    en_t.θ_ρ = virtual_pottemp(ts)
end;

# Specify where in `diffusive::Vars` to store the computed gradient from
function compute_gradient_flux!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT, N}
    # println("-------- compute_gradient_flux!")
    # Aliases:
    gm = state
    gm_d = diffusive
    up_d = diffusive.turbconv.updraft
    up_∇t = ∇transform.turbconv.updraft
    en_d = diffusive.turbconv.environment
    en_∇t = ∇transform.turbconv.environment

    for i in 1:N
        up_d[i].∇u = up_∇t[i].u
    end

    ρinv = 1 / gm.ρ
    # negative signs here as we have a '-' sign in BL form leading to + K∂ϕ/∂z on the RHS
    # first moment grid mean comming from enviroment gradients only
    en_d.∇e = en_∇t.e
    en_d.∇e_int = en_∇t.e_int
    en_d.∇q_tot = en_∇t.q_tot
    en_d.∇u = en_∇t.u
    # second moment env cov
    en_d.∇tke = en_∇t.tke
    en_d.∇e_int_cv = en_∇t.e_int_cv
    en_d.∇q_tot_cv = en_∇t.q_tot_cv
    en_d.∇e_int_q_tot_cv = en_∇t.e_int_q_tot_cv

    en_d.∇θ_ρ = en_∇t.θ_ρ

end;

# We have no sources, nor non-diffusive fluxes.

function turbconv_source!(
    m::AtmosModel{FT},
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
) where {FT}
    # println("-------- turbconv_source!")

    turbconv = m.turbconv
    N = n_updrafts(m.turbconv)
    # Aliases:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    gm_s = source
    en_s = source.turbconv.environment
    up_s = source.turbconv.updraft
    en_d = diffusive.turbconv.environment


    # grid mean sources - I think that large scale subsidence in
    #            doubly periodic domains should be applied here
    # updraft sources

    # YAIR  - these need to be defined as vectors length N - check with Charlie
    εt = MArray{Tuple{N}, FT}(zeros(FT, N))
    ε = MArray{Tuple{N}, FT}(zeros(FT, N))
    δ = MArray{Tuple{N}, FT}(zeros(FT, N))
    # should be conditioned on updraft_area > minval

    # get environment values for e, q_tot , u[3]
    _grav::FT = grav(m.param_set)
    ρinv = 1 / gm.ρ
    a_env = 1 - sum([up[i].ρa for i in 1:N]) * ρinv
    w_env = (gm.ρu[3] - sum([up[i].ρau[3] for i in 1:N])) * ρinv
    e_env = (gm.ρe - sum([up[i].ρae for i in 1:N])) * ρinv
    q_tot_env = (gm.moisture.ρq_tot - sum([up[i].ρaq_tot for i in 1:N])) * ρinv

    u_env = (gm.ρu .- sum([up[i].ρau for i in 1:N])) / (gm.ρ * a_env)

    en_ρe = (gm.ρe - sum([up[j].ρae for j in 1:N])) / a_env
    en_ρu = (gm.ρu .- sum([up[j].ρae for j in 1:N])) / a_env
    e_pot = gravitational_potential(m, aux)
    e_int_env = internal_energy(gm.ρ, en_ρe, en_ρu, e_pot)
    e_int_gm = internal_energy(gm.ρ, gm.ρe, gm.ρu, e_pot)

    for i in 1:N
        # upd vars
        w_i = up[i].ρau[3] / up[i].ρa
        e_int_upd = internal_energy(
            gm.ρ,
            up[i].ρae / up[i].ρa,
            up[i].ρau / up[i].ρa,
            e_pot,
        )

        # first moment sources
        # ε[i], δ[i], εt[i] = entr_detr(m, turbconv.entr_detr, state, aux, t, i)
        # dpdz, dpdz_tke_i  = perturbation_pressure(m, turbconv.pressure, state, diffusive, aux, t, direction, i)
        ε[i] = FT(0)
        δ[i] = FT(0)
        εt[i] = FT(0)
        dpdz = FT(0)
        dpdz_tke_i = FT(0)

        # entrainment and detrainment
        up_s[i].ρa += up[i].ρa * w_i * (ε[i] - δ[i])
        up_s[i].ρau +=
            up[i].ρa *
            w_i *
            ((ε[i] + εt[i]) * u_env - (δ[i] + εt[i]) * up_s[i].ρau)
        up_s[i].ρae +=
            up[i].ρa *
            w_i *
            ((ε[i] + εt[i]) * e_env - (δ[i] + εt[i]) * up_s[i].ρae)
        up_s[i].ρaq_tot +=
            up[i].ρa *
            w_i *
            ((ε[i] + εt[i]) * q_tot_env - (δ[i] + εt[i]) * up_s[i].ρaq_tot)

        # perturbation pressure in w equation
        up_s[i].ρau += SVector(0, 0, up[i].ρa * dpdz)

        # microphysics sources should be applied here

        ## environment second moments:

        # covariances entrinament sources from the i'th updraft
        # need to compute e_int in updraft and gridmean for entrainment
        # -- if ϕ'ψ' is tke and ϕ,ψ are both w than a factor 0.5 appears in the εt and δ terms
        # Covar_Source      +=  ρaw⋅δ⋅(ϕ_up-ϕ_en) ⋅ (ψ_up-ψ_en)
        #                     + ρaw⋅εt⋅[(ϕ_up-⟨ϕ⟩)⋅(ψ_up-ψ_en)
        #                     + (ϕ_up-⟨ϕ⟩)⋅(ψ_up-ψ_en)] - ρaw⋅ε⋅ϕ'ψ'

        en_s.ρatke += (
            up[i].ρau[3] *
            δ[i] *
            (up[i].ρau[3] / up[i].ρa - w_env) *
            (up[i].ρau[3] / up[i].ρa - w_env) *
            FT(0.5) +
            up[i].ρau[3] *
            εt[i] *
            w_env *
            (up[i].ρau[3] / up[i].ρa - gm.ρu[i] * ρinv) -
            up[i].ρau[3] * ε[i] * en.ρatke
        )

        # en_s.ρae_cv += (
        #     up[i].ρau[3] *
        #     δ[i] *
        #     (e_int_upd - e_int_env) *
        #     (e_int_upd - e_int_env) +
        #     up[i].ρau[3] * εt[i] * e_int_env * (e_int_upd - e_int_gm) * 2 -
        #     up[i].ρau[3] * ε[i] * en.ρae_int_cv
        # )

        en_s.ρaq_tot_cv += (
            up[i].ρau[3] *
            δ[i] *
            (up[i].ρaq_tot / up[i].ρa - q_tot_env) *
            (up[i].ρaq_tot / up[i].ρa - q_tot_env) +
            up[i].ρau[3] *
            εt[i] *
            q_tot_env *
            (up[i].ρaq_tot / up[i].ρa - gm.moisture.ρq_tot * ρinv) *
            2 - up[i].ρau[3] * ε[i] * en.ρaq_tot_cv
        )

        # en_s.ρae_q_tot_cv += (
        #     up[i].ρau[3] *
        #     δ[i] *
        #     (e_int_upd - e_int_env) *
        #     (up[i].ρaq_tot / up[i].ρa - q_tot_env) +
        #     up[i].ρau[3] * εt[i] * q_tot_env * (e_int_upd - e_int_gm) +
        #     up[i].ρau[3] *
        #     εt[i] *
        #     e_int_env *
        #     (up[i].ρaq_tot / up[i].ρa - gm.moisture.ρq_tot * ρinv) -
        #     up[i].ρau[3] * ε[i] * en.ρae_q_tot_cv
        # )

        # pressure tke source from the i'th updraft
        en_s.ρatke += up[i].ρa * dpdz_tke_i
    end
    # mix_len    = mixing_length(m, turbconv.mix_len, state, diffusive, aux, t, δ, εt)
    mix_len = FT(0)
    en_tke = en.ρatke * ρinv / a_env
    K_eddy = turbconv.mix_len.c_k * mix_len * sqrt(en_tke)
    Shear = en_d.∇u[1, 3] .^ 2 + en_d.∇u[2, 3] .^ 2 + en_d.∇u[3, 3] .^ 2 # consider scalar product of two vectors

    # second moment production from mean gradients (+ sign here as we have + S in BL form)
    #                            production from mean gradient           - Dissipation
    en_s.ρatke += gm.ρ * a_env * K_eddy * Shear
    -turbconv.mix_len.c_m * sqrt(en_tke) / mix_len * en.ρatke
    en_s.ρae_int_cv += gm.ρ * a_env * K_eddy * en_d.∇e_int[3] * en_d.∇e_int[3]
    -turbconv.mix_len.c_m * sqrt(en_tke) / mix_len * en.ρae_int_cv
    en_s.ρaq_tot_cv += gm.ρ * a_env * K_eddy * en_d.∇q_tot[3] * en_d.∇q_tot[3]
    -turbconv.mix_len.c_m * sqrt(en_tke) / mix_len * en.ρaq_tot_cv
    en_s.ρae_int_q_tot_cv +=
        gm.ρ * a_env * K_eddy * en_d.∇e_int[3] * en_d.∇q_tot[3]
    -turbconv.mix_len.c_m * sqrt(en_tke) / mix_len * en.ρae_int_q_tot_cv
    # covariance microphysics sources should be applied here
end;

# in the EDMF first order (advective) fluxes exist only in the grid mean (if <w> is nonzero) and the uprdafts
function flux_first_order!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT, N}
    # println("-------- flux_first_order!")

    # Aliases:
    gm = state
    up = state.turbconv.updraft
    up_f = flux.turbconv.updraft

    # positive sign here as we have a '-' sign in BL form leading to - ∂ρwϕ/∂z on the RHS
    # updrafts
    ρinv = 1 / gm.ρ
    for i in 1:N
        up_f[i].ρa = up[i].ρau
        u = up[i].ρau / up[i].ρa
        up_f[i].ρau = up[i].ρau * u'
        up_f[i].ρae = u * up[i].ρae
        up_f[i].ρaq_tot = u * up[i].ρaq_tot
    end

end;

# in the EDMF second order (diffusive) fluxes exist only in the grid mean and the enviroment
function flux_second_order!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
) where {FT, N}
    # println("-------- flux_second_order!")
    # Aliases:
    gm = state
    up = state.turbconv.updraft
    en = state.turbconv.environment
    gm_f = flux
    up_f = flux.turbconv.updraft
    en_f = flux.turbconv.environment
    en_d = diffusive.turbconv.environment

    ρinv = FT(1) / gm.ρ
    en_ρa = gm.ρ - sum([up[i].ρa for i in 1:N])

    εt = MArray{Tuple{N}, FT}(zeros(FT, N))
    ε = MArray{Tuple{N}, FT}(zeros(FT, N))
    δ = MArray{Tuple{N}, FT}(zeros(FT, N))
    for i in 1:N
        # ε[i], δ[i], εt[i] = entr_detr(m, m.turbconv.entr_detr, state, aux, t, i)
        ε[i] = FT(0)
        δ[i] = FT(0)
        εt[i] = FT(0)
    end
    # mix_len = mixing_length(m, turbconv.mix_len, state, diffusive, aux, t, δ, εt)
    mix_len = FT(0)
    ρa_env = (gm.ρ - sum([up[j].ρa for j in 1:N]))
    K_eddy = m.turbconv.mix_len.c_k * mix_len * sqrt(abs(en.ρatke / ρa_env)) # YAIR

    ## we are adding the massflux term here as it is part of the total flux:
    #total flux(ϕ) =   diffusive_flux(ϕ)  +   massflux(ϕ)
    #   ⟨w ⃰ ϕ ⃰ ⟩   = - a_0 K_eddy⋅∂ϕ/∂z + ∑ a_i(w_i-⟨w⟩)(ϕ_i-⟨ϕ⟩)
    massflux_e = sum([
        up[i].ρa *
        ρinv *
        (gm.ρe * ρinv - up[i].ρae / up[i].ρa) *
        (gm.ρu[3] * ρinv - up[i].ρau[3] / up[i].ρa) for i in 1:N
    ])
    massflux_q_tot = sum([
        up[i].ρa *
        ρinv *
        (gm.moisture.ρq_tot * ρinv - up[i].ρaq_tot / up[i].ρa) *
        (gm.ρu[3] * ρinv - up[i].ρau[3] / up[i].ρa) for i in 1:N
    ])
    massflux_u = sum([
        up[i].ρa *
        ρinv *
        (gm.ρu * ρinv .- up[i].ρau / up[i].ρa) *
        (gm.ρu[3] * ρinv - up[i].ρau[3] / up[i].ρa) for i in 1:N
    ])

    # update grid mean flux_second_order
    gm_f.ρe += -ρa_env * K_eddy * en_d.∇e[3] + massflux_e
    gm_f.moisture.ρq_tot += -ρa_env * K_eddy * en_d.∇q_tot[3] + massflux_q_tot
    # override ---------------- adding massflux below u breaks
    # gm_f.ρu     += -ρa_env.*K_eddy.*en_d.∇u
    gm_f.ρu += SMatrix{3, 3, FT, 9}(
        0,
        0,
        0,
        0,
        0,
        0,
        -ρa_env .* K_eddy .* en_d.∇u[1, 3] + massflux_u[1],
        -ρa_env .* K_eddy .* en_d.∇u[2, 3] + massflux_u[2],
        -ρa_env .* K_eddy .* en_d.∇u[3, 3] + massflux_u[3],
    )

    # env second momment flux_second_order
    en_f.ρatke += -ρa_env * K_eddy * en_d.∇tke[3]
    en_f.ρae_int_cv += -ρa_env * K_eddy * en_d.∇e_int_cv[3]
    en_f.ρaq_tot_cv += -ρa_env * K_eddy * en_d.∇q_tot_cv[3]
    en_f.ρae_int_q_tot_cv += -ρa_env * K_eddy * en_d.∇e_int_q_tot_cv[3]
end;

# ### Boundary conditions

# Second-order terms in our equations, ``∇⋅(G)`` where ``G = α∇ρcT``, are
# internally reformulated to first-order unknowns.
# Boundary conditions must be specified for all unknowns, both first-order and
# second-order unknowns which have been reformulated.

# The boundary conditions for `ρaq_tot` (first order unknown)
function turbconv_boundary_state!(
    nf,
    bc::EDMFBCs,
    m::AtmosModel{FT},
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
) where {FT}
    # println("-------- turbconv_boundary_state!")
    turbconv = m.turbconv
    N = n_updrafts(turbconv)
    up = state⁺.turbconv.updraft
    en = state⁺.turbconv.environment
    gm = state⁺
    if bctype == 1 # bottom
        # YAIR - questions which state should I use here , state⁺ or state⁻  for computation of surface processes
        upd_a_surf, upd_e_surf, upd_q_tot_surf =
            compute_updraft_surface_BC(turbconv.surface, turbconv, m, gm)
        upd_a_surf = FT(0.1)
        upd_e_surf = FT(210000)
        upd_q_tot_surf = FT(0.0016)
        for i in 1:N
            up[i].ρau = SVector(0, 0, 0)
            up[i].ρa = upd_a_surf
            up[i].ρae = upd_e_surf
            up[i].ρaq_tot = upd_q_tot_surf
        end

        tke, e_cv, q_tot_cv, e_q_tot_cv =
            env_surface_covariances(turbconv.surface, turbconv, m, gm)
        area_en = 1 - sum([up[i].ρa for i in 1:N]) / gm.ρ
        tke = FT(0)
        e_int_cv = FT(0)
        q_tot_cv = FT(0)
        e_int_q_tot_cv = FT(0)
        en.ρatke = gm.ρ * area_en * tke
        en.ρae_int_cv = gm.ρ * area_en * e_cv
        en.ρaq_tot_cv = gm.ρ * area_en * q_tot_cv
        en.ρae_int_q_tot_cv = gm.ρ * area_en * e_q_tot_cv

    elseif bctype == 2 # top
        ρinv = 1 / gm.ρ
        for i in 1:N
            up[i].ρau = SVector(0, 0, 0)
            up[i].ρa = FT(0)
            up[i].ρae = FT(0)
            up[i].ρaq_tot = FT(0)
        end
        en.ρatke = FT(0)
        en.ρae_int_cv = FT(0)
        en.ρaq_tot_cv = FT(0)
        en.ρae_int_q_tot_cv = FT(0)
    end
end;

# The boundary conditions for `ρcT` are specified here for second-order
# unknowns
function turbconv_normal_boundary_flux_second_order!(
    nf,
    bc::EDMFBCs,
    m::AtmosModel{FT},
    fluxᵀn::Vars,
    n⁻,
    state⁻::Vars,
    diff⁻::Vars,
    hyperdiff⁻::Vars,
    aux⁻::Vars,
    state⁺::Vars,
    diff⁺::Vars,
    hyperdiff⁺::Vars,
    aux⁺::Vars,
    bctype,
    t,
    _...,
) where {FT}
    # println("-------- turbconv_normal_boundary_flux_second_order!")
    turbconv = m.turbconv
    N = n_updrafts(turbconv)
    # gm = state⁺
    # up = state⁺.turbconv.updraft
    # en = state⁺.turbconv.environment
    # gm_d = diff⁺
    # up_d = diff⁺.turbconv.updraft
    # en_d = diff⁺.turbconv.environment
    # # Charlie is the use of state⁺ here consistent for gm.ρ, up[i].ρa ?
    # if bctype == 1 # bottom
    #     area_en  = 1 - sum([up[i].ρa for i in 1:N])/gm.ρ
    #     # YAIR - I need to pass the SurfaceModel into BC and into env_surface_covariances
    #     # tke, e_cv ,q_tot_cv ,e_q_tot_cv = env_surface_covariances(turbconv.surface, turbconv, m, gm)
    #     tke = FT(0)
    #     e_int_cv = FT(0)
    #     q_tot_cv = FT(0)
    #     e_int_q_tot_cv = FT(0)
    #     en_d.∇tke = SVector(0,0,gm.ρ * area_en * tke)
    #     en_d.∇e_int_cv = SVector(0,0,gm.ρ * area_en * e_cv)
    #     en_d.∇q_tot_cv = SVector(0,0,gm.ρ * area_en * q_tot_cv)
    #     en_d.∇e_int_q_tot_cv = SVector(0,0,gm.ρ * area_en * e_q_tot_cv)
    # elseif bctype == 2 # top
    #     # for now zero flux at the top
    #     en_d.∇tke = -n⁻ * FT(0)
    #     en_d.∇e_int_cv = -n⁻ * FT(0)
    #     en_d.∇q_tot_cv = -n⁻ * FT(0)
    #     en_d.∇e_int_q_tot_cv = -n⁻ * FT(0)
    # end
end;
