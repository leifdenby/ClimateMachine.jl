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

# For debugging
debug_kernels = false
kernel_calls = Dict([
    :init_state_conservative! => false,
    :init_aux_turbconv! => false,
    :update_auxiliary_state! => false,
    :flux_first_order! => false,
    :flux_second_order! => false,
    :turbconv_boundary_state! => false,
    :turbconv_normal_boundary_flux_second_order! => false,
    :compute_gradient_flux! => false,
])

using Printf
import ClimateMachine.TurbulenceConvection:
    vars_state_auxiliary,
    vars_state_conservative,
    vars_state_gradient,
    vars_state_gradient_flux,
    init_aux_turbconv!,
    update_auxiliary_state!,
    flux_first_order!,
    flux_second_order!,
    turbconv_sources,
    turbconv_boundary_state!,
    turbconv_normal_boundary_flux_second_order!,
    compute_gradient_argument!,
    compute_gradient_flux!
import ..Thermodynamics: air_pressure, air_density


include(joinpath("helper_funcs", "nondimensional_exchange_functions.jl"))
include(joinpath("helper_funcs", "lamb_smooth_minimum.jl"))
include(joinpath("helper_funcs", "subdomain_statistics.jl"))
include(joinpath("helper_funcs", "diagnose_environment.jl"))
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
    @vars(ρa::FT, ρau::SVector{3, FT}, ρaθ_liq::FT, ρaq_tot::FT,)
end

function vars_state_conservative(::Environment, FT)
    @vars(ρatke::FT, ρaθ_liq_cv::FT, ρaq_tot_cv::FT, ρaθ_liq_q_tot_cv::FT,)
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
        θ_liq::FT,
        q_tot::FT,
        u::SVector{3, FT},
        tke::FT,
        θ_liq_cv::FT,
        q_tot_cv::FT,
        θ_liq_q_tot_cv::FT,
        θ_vap::FT,
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
        ∇θ_liq::SVector{3, FT},
        ∇q_tot::SVector{3, FT},
        ∇u::SMatrix{3, 3, FT, 9},
        ∇tke::SVector{3, FT},
        ∇θ_liq_cv::SVector{3, FT},
        ∇q_tot_cv::SVector{3, FT},
        ∇θ_liq_q_tot_cv::SVector{3, FT},
        ∇θ_vap::SVector{3, FT},# used in a diagnostic equation for the mixing length
    )

end

function vars_state_gradient_flux(m::EDMF, FT)
    @vars(
        environment::vars_state_gradient_flux(m.environment, FT),
        updraft::vars_state_gradient_flux(m.updraft, FT)
    )
end

function init_aux_turbconv!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    aux::Vars,
    geom::LocalGeometry,
) where {FT, N}
    kernel_calls[:init_aux_turbconv!] = true
    debug_kernels && println("Calling init_aux_turbconv!")
    debug_kernels && @show kernel_calls

    # # Aliases:
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft

    en_a.cld_frac = eps(FT)

    for i in 1:N
        up_a[i].buoyancy = eps(FT)
        up_a[i].updraft_top = 500 #eps(FT)
    end
    en_a.cld_frac = eps(FT)
end;

# - this method is only called at `t=0`
function init_state_conservative!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
) where {FT, N}
    kernel_calls[:init_state_conservative!] = true
    debug_kernels && println("Calling init_state_conservative!")
    debug_kernels && @show kernel_calls

    # # Aliases:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft

    # GCM setting - Initialize the grid mean profiles of prognostic variables (ρ,e_int,q_tot,u,v,w)
    z = altitude(m, aux)

    # SCM setting - need to have separate cases coded and called from a folder - see what LES does
    # a moist_thermo state is used here to convert the input θ,q_tot to e_int, q_tot profile
    @info(state.ρ , state.ρe, state.mositure.ρq_tot)
    ts = thermo_state(m, state, aux)
    θ_liq = liquid_ice_pottemp(ts)
    θ_liq = FT(280)

    a_up = FT(0.1)
    for i in 1:N
        up[i].ρa = gm.ρ * a_up
        up[i].ρau = gm.ρu * a_up
        up[i].ρaθ_liq = gm.ρ * a_up * θ_liq
        up[i].ρaq_tot = gm.moisture.ρq_tot * a_up
    end

    # initialize environment covariance with zero for now
    en.ρaθ_liq_cv = eps(FT)
    en.ρaq_tot_cv = eps(FT)
    en.ρaθ_liq_q_tot_cv = eps(FT)
    return nothing
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
    kernel_calls[:update_auxiliary_state!] = true
    debug_kernels && println("Calling update_auxiliary_state!")
    debug_kernels && @show kernel_calls

    turbconv = m.turbconv

    # ----------
    # turbconv.updraft[i].updraft_top should be
    # defined after this method call, but it's
    # not complete yet...
    # compute_updraft_top!(dg, m, turbconv, Q, t, elems)
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
    kernel_calls[:edmf_stack_nodal_update_aux!] = true
    debug_kernels && println("Calling edmf_stack_nodal_update_aux!")
    debug_kernels && @show kernel_calls

    N = n_updrafts(m.turbconv)

    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft

    #  -------------  Compute buoyancies of subdomains
    ts = thermo_state(m, state, aux)
    gm_p = air_pressure(ts)
    ρinv = 1 / gm.ρ
    _grav::FT = grav(m.param_set)

    for i in 1:N
        # T_i = saturation_adjustment_q_tot_θ_liq_ice_given_pressure(m.param_set, up[i].ρaθ_liq*ρinv, gm_p, up[i].ρq_tot*ρinv, PhaseEquil, Int64(10), FT(1e-2))
        # ρ_i = air_density(param_set, T_i, gm_p, PhasePartition(up[i].ρaq_tot*ρinv))
        # ts = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, up[i].ρaθ_liq/up[i].ρa, gm_p, up[i].ρaq_tot/up[i].ρa)
        # ρ_i = air_density(ts)
        # override
        ρ_i = gm.ρ
        up_a[i].buoyancy = -_grav * (ρ_i - aux.ref_state.ρ) * ρinv
    end
    # compute the buoyancy of the environment
    en_θ_liq = environment_θ_liq(m, state, aux, N)
    en_q_tot = environment_q_tot(state, aux, N)
    # ts = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, en_θ_liq, gm_p, en_q_tot)
    # en_ρ = air_density(ts)
    en_ρ = gm.ρ
    b_env = -_grav * (en_ρ - aux.ref_state.ρ) * ρinv
    en_area = environment_area(state, aux, N)
    b_gm = en_area * b_env + sum([up_a[i].buoyancy*up[i].ρa*ρinv for i in 1:N])
    # loop over all updrafts and remove the gm_b
    for i in 1:N
        up_a[i].buoyancy -= b_gm
    end
end;

enforce_unit_bounds(x) = clamp(x, 1e-3, 1-1e-3)
# enforce_unit_bounds(x) = x
enforce_positivity(x) = max(x, 0)
# enforce_positivity(x) = x

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
    kernel_calls[:compute_gradient_argument!] = true
    debug_kernels && println("Calling compute_gradient_argument!")
    debug_kernels && @show kernel_calls

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
    en_area   = environment_area(state,aux,N)

    # populate gradient arguments
    en_t.θ_liq = environment_θ_liq(m, state ,aux ,N)
    en_t.q_tot = environment_q_tot(state ,aux ,N)
    en_t.u     = environment_u(state ,aux ,N)

    en_t.tke            = enforce_positivity(en.ρatke / (en_area * gm.ρ))
    en_t.θ_liq_cv       = en.ρaθ_liq_cv / (en_area * gm.ρ)
    en_t.q_tot_cv       = en.ρaq_tot_cv / (en_area * gm.ρ)
    en_t.θ_liq_q_tot_cv = en.ρaθ_liq_q_tot_cv / (en_area * gm.ρ)

    en_t.θ_vap = virtual_pottemp(ts)
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
    kernel_calls[:compute_gradient_flux!] = true
    debug_kernels && println("Calling compute_gradient_flux!")
    debug_kernels && @show kernel_calls

    # # Aliases:
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
    en_d.∇θ_liq = en_∇t.θ_liq
    en_d.∇q_tot = en_∇t.q_tot
    en_d.∇u = en_∇t.u
    # second moment env cov
    en_d.∇tke = en_∇t.tke
    en_d.∇θ_liq_cv = en_∇t.θ_liq_cv
    en_d.∇q_tot_cv = en_∇t.q_tot_cv
    en_d.∇θ_liq_q_tot_cv = en_∇t.θ_liq_q_tot_cv

    en_d.∇_vap = en_∇t.θ_vap
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
    kernel_calls[:turbconv_source!] = true
    debug_kernels && println("Calling turbconv_source!")
    debug_kernels && @show kernel_calls

    # turbconv = m.turbconv
    N = n_updrafts(m.turbconv)
    # Aliases:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    gm_s = source
    en_s = source.turbconv.environment
    up_s = source.turbconv.updraft
    en_d = diffusive.turbconv.environment
    ρatke_env = enforce_positivity(en.ρatke)


    # grid mean sources - I think that large scale subsidence in
    #            doubly periodic domains should be applied here
    εt = MArray{Tuple{N}, FT}(zeros(FT, N))
    ε = MArray{Tuple{N}, FT}(zeros(FT, N))
    δ = MArray{Tuple{N}, FT}(zeros(FT, N))

    # get environment values for e, q_tot , u[3]
    _grav::FT = grav(m.param_set)
    ρinv = 1 / gm.ρ
    en_a     = environment_area(state, aux, N)
    en_w     = environment_w(state, aux, N)
    en_θ_liq = environment_θ_liq(m, state, aux, N)
    en_q_tot = environment_q_tot(state, aux, N)
    en_u     = environment_u(state, aux, N)
    ts = thermo_state(m, state, aux)
    # gm_θ_liq = liquid_ice_pottemp(ts)
    gm_θ_liq = FT(280)

    # en_ρu     = gm.ρ.* en_u
    # en_ρθ_liq = gm.ρ * en_ρθ_liq
    # en_ρq_tot = gm.ρ * en_ρq_tot

    for i in 1:N
        # upd vars
        w_i = up[i].ρau[3] / up[i].ρa
        ρa_i = enforce_unit_bounds(up[i].ρa)

        # first moment sources
        # ε_dyn[i] ,δ_dyn[i], ε_trb[i] = entr_detr(m, turbconv.entr_detr, state, aux, t, i)
        # dpdz, dpdz_tke_i  = perturbation_pressure(m, turbconv.pressure, state, diffusive, aux, t, direction, i)
        ε_dyn[i] = FT(0)
        δ_dyn[i] = FT(0)
        ε_trb[i] = FT(0)
        dpdz = FT(0)
        dpdz_tke_i = FT(0)

        # entrainment and detrainment
        up_s[i].ρa  += up[i].ρau[3] * (ε_dyn[i] - δ_dyn[i])
        up_s[i].ρau += up[i].ρau[3] * # YAIR is this correct?
            ((ε_dyn[i] + ε_trb[i]) * en_u .- (δ_dyn[i] + ε_trb[i]) * up[i].ρau/ρa_i)
        up_s[i].ρaθ_liq += up[i].ρau[3] *
            ((ε_dyn[i] + ε_trb[i]) * en_θ_liq - (δ_dyn[i] + ε_trb[i]) * up[i].ρaθ_liq/ρa_i)
        up_s[i].ρaq_tot += up[i].ρau[3] *
            ((ε_dyn[i] + ε_trb[i]) * en_q_tot - (δ_dyn[i] + ε_trb[i]) * up[i].ρaq_tot/ρa_i)

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
            δ_dyn[i] *
            (up[i].ρau[3]/ρa_i - en_w) *
            (up[i].ρau[3]/ρa_i - en_w) *
            FT(0.5) +
            up[i].ρau[3] *
            ε_trb[i] *
            en_w *
            (up[i].ρau[3]/ρa_i - gm.ρu[3]*ρinv) -
            up[i].ρau[3]*ε_dyn[i] * ρatke_env
        )

        en_s.ρaθ_liq_cv += (
            up[i].ρau[3] *
            δ_dyn[i] *
            (up[i].ρaθ_liq/ρa_i - en_θ_liq) *
            (up[i].ρaθ_liq/ρa_i - en_θ_liq) +
            2 * up[i].ρau[3] * ε_trb[i] * en_θ_liq * (up[i].ρaθ_liq/ρa_i - gm_θ_liq) -
            up[i].ρau[3] * ε_dyn[i]  * en.ρaθ_liq_cv
        )

        en_s.ρaq_tot_cv += (
            up[i].ρau[3] *
            δ_dyn[i] *
            (up[i].ρaq_tot/ρa_i - en_q_tot) *
            (up[i].ρaq_tot/ρa_i - en_q_tot) +
            2 * up[i].ρau[3] * ε_trb[i] * en_q_tot *
            (up[i].ρaq_tot / ρa_i - gm.moisture.ρq_tot * ρinv) -
            up[i].ρau[3] * ε_dyn[i] * en.ρaq_tot_cv
        )

        en_s.ρaθ_liq_q_tot_cv += (
            up[i].ρau[3] *
            δ_dyn[i] *
            (up[i].ρaθ_liq/ρa_i - en_θ_liq) *
            (up[i].ρaq_tot/ρa_i - en_q_tot) +
            up[i].ρau[3] * ε_trb[i] * en_q_tot     *
            (up[i].ρaθ_liq/ρa_i - gm_θ_liq) +
            up[i].ρau[3] * ε_trb[i] * en_θ_liq     *
            (up[i].ρaq_tot/ρa_i - gm.moisture.ρq_tot * ρinv) -
            up[i].ρau[3] * ε_dyn[i]  * en.ρaθ_liq_q_tot_cv
        )

        # pressure tke source from the i'th updraft
        en_s.ρatke += up[i].ρa * dpdz_tke_i
    end
    # l_mix    = mixing_length(m, turbconv.mix_len, state, diffusive, aux, t, δ, εt)
    l_mix = FT(0)
    # en_tke = ρatke_env * ρinv / en_a
    K_eddy = m.turbconv.mix_len.c_k * l_mix * sqrt(ρatke_env)
    Shear = en_d.∇u[1, 3] .^ 2 + en_d.∇u[2, 3] .^ 2 + en_d.∇u[3, 3] .^ 2 # consider scalar product of two vectors

    # second moment production from mean gradients (+ sign here as we have + S in BL form)
    #                            production from mean gradient           - Dissipation
    en_s.ρatke            += gm.ρ * en_a * K_eddy * Shear
    -m.turbconv.mix_len.c_m * sqrt(ρatke_env) / l_mix * ρatke_env
    en_s.ρaθ_liq_cv       += gm.ρ * en_a * K_eddy * en_d.∇θ_liq[3] * en_d.∇θ_liq[3]
    -m.turbconv.mix_len.c_m * sqrt(ρatke_env) / l_mix * en.ρaθ_liq_cv
    en_s.ρaq_tot_cv       += gm.ρ * en_a * K_eddy * en_d.∇q_tot[3] * en_d.∇q_tot[3]
    -m.turbconv.mix_len.c_m * sqrt(ρatke_env) / l_mix * en.ρaq_tot_cv
    en_s.ρaθ_liq_q_tot_cv +=  gm.ρ * en_a * K_eddy * en_d.∇θ_liq[3] * en_d.∇q_tot[3]
    -m.turbconv.mix_len.c_m * sqrt(ρatke_env) / l_mix * en.ρaθ_liq_q_tot_cv
#     # covariance microphysics sources should be applied here
end;

# # in the EDMF first order (advective) fluxes exist only in the grid mean (if <w> is nonzero) and the uprdafts
function flux_first_order!(
    turbconv::EDMF{FT, N},
    m::AtmosModel{FT},
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT, N}
    kernel_calls[:flux_first_order!] = true
    debug_kernels && println("Calling flux_first_order!")
    debug_kernels && @show kernel_calls

    # # Aliases:
    gm = state
    up = state.turbconv.updraft
    up_f = flux.turbconv.updraft

    # positive sign here as we have a '-' sign in BL form leading to - ∂ρwϕ/∂z on the RHS
    # updrafts
    ρinv = 1 / gm.ρ
    for i in 1:N
        ρa_i = enforce_unit_bounds(up[i].ρa)
        up_f[i].ρa = up[i].ρau
        u = up[i].ρau / ρa_i
        up_f[i].ρau = up[i].ρau * u'
        up_f[i].ρaθ_liq = u * up[i].ρaθ_liq
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
    kernel_calls[:flux_second_order!] = true
    debug_kernels && println("Calling flux_second_order!")
    debug_kernels && @show kernel_calls

    # Aliases:
    gm = state
    up = state.turbconv.updraft
    en = state.turbconv.environment
    gm_f = flux
    up_f = flux.turbconv.updraft
    en_f = flux.turbconv.environment
    en_d = diffusive.turbconv.environment
    ρatke_env = enforce_positivity(en.ρatke)
    ρinv = FT(1) / gm.ρ

    εt = MArray{Tuple{N}, FT}(zeros(FT, N))
    ε = MArray{Tuple{N}, FT}(zeros(FT, N))
    δ = MArray{Tuple{N}, FT}(zeros(FT, N))
    for i in 1:N
        # ε[i], δ[i], εt[i] = entr_detr(m, m.turbconv.entr_detr, state, aux, t, i)
        ε[i] = FT(0)
        δ[i] = FT(0)
        εt[i] = FT(0)
    end
    # l_mix = mixing_length(m, turbconv.mix_len, state, diffusive, aux, t, δ, εt)
    l_mix = FT(0)
    en_area = environment_area(state, aux, N)
    K_eddy = m.turbconv.mix_len.c_k * l_mix * sqrt(abs(ρatke_env / en_area * ρinv)) # YAIR

    ## we are adding the massflux term here as it is part of the total flux:
    #total flux(ϕ) =   diffusive_flux(ϕ)  +   massflux(ϕ)
    #   ⟨w ⃰ ϕ ⃰ ⟩   = - a_0 K_eddy⋅∂ϕ/∂z + ∑ a_i(w_i-⟨w⟩)(ϕ_i-⟨ϕ⟩)

    massflux_e = FT(0)
    ts = thermo_state(m, state, aux)
    gm_p = air_pressure(ts)
    gm_p = FT(100000)
    for i in 1:N
        # ts = LiquidIcePotTempSHumEquil_given_pressure(m.param_set, up[i].ρaθ_liq/up[i].ρa, gm_p, up[i].ρaq_tot/up[i].ρa)
        # e_kin = FT(1 // 2) * (up[i].u[1]^2 + up[i].u[2]^2 + up[i].u[3]^2)
        # up_e = total_energy(e_kin, _grav * z, ts)
        ρa_i = enforce_unit_bounds(up[i].ρa)
        up_e = gm.ρe * ρinv
        massflux_e += up[i].ρa *
            ρinv * (gm.ρe * ρinv - up_e) *
            (gm.ρu[3] * ρinv - up[i].ρau[3] / ρa_i)
    end

    massflux_q_tot = sum([
        up[i].ρa *
        ρinv *
        (gm.moisture.ρq_tot * ρinv - up[i].ρaq_tot / up[i].ρa) *
        (gm.ρu[3] * ρinv - up[i].ρau[3] / enforce_unit_bounds(up[i].ρa)) for i in 1:N
    ])
    massflux_u = sum([
        up[i].ρa *
        ρinv *
        (gm.ρu * ρinv .- up[i].ρau / up[i].ρa) *
        (gm.ρu[3] * ρinv - up[i].ρau[3] / enforce_unit_bounds(up[i].ρa)) for i in 1:N
    ])

    # update grid mean flux_second_order

    # This was commented for debugging purposes:
    ∂e∂θl = FT(1)

    gm_f.ρe              += - gm.ρ*en_area * K_eddy * en_d.∇θ_liq[3]*∂e∂θl + massflux_e
    gm_f.moisture.ρq_tot += - gm.ρ*en_area * K_eddy * en_d.∇q_tot[3] + massflux_q_tot
    gm_f.ρu = gm_f.ρu .+ SMatrix{3, 3, FT, 9}(
        0,
        0,
        0,
        0,
        0,
        0,
        -gm.ρ*en_area .* K_eddy .* en_d.∇u[1, 3] + massflux_u[1],
        -gm.ρ*en_area .* K_eddy .* en_d.∇u[2, 3] + massflux_u[2],
        -gm.ρ*en_area .* K_eddy .* en_d.∇u[3, 3] + massflux_u[3],
    )

    # env second momment flux_second_order
    en_f.ρatke            += -gm.ρ*en_area * K_eddy * en_d.∇tke[3]
    en_f.ρaθ_liq_cv       += -gm.ρ*en_area * K_eddy * en_d.∇θ_liq_cv[3]
    en_f.ρaq_tot_cv       += -gm.ρ*en_area * K_eddy * en_d.∇q_tot_cv[3]
    en_f.ρaθ_liq_q_tot_cv += -gm.ρ*en_area * K_eddy * en_d.∇θ_liq_q_tot_cv[3]
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
    kernel_calls[:turbconv_boundary_state!] = true
    debug_kernels && println("Calling turbconv_boundary_state!")
    debug_kernels && @show kernel_calls

    turbconv = m.turbconv
    N = n_updrafts(turbconv)
    up = state⁺.turbconv.updraft
    en = state⁺.turbconv.environment
    gm = state⁺
    gm_a = aux⁺
    if bctype == 1 # bottom
        # YAIR - questions which state should I use here , state⁺ or state⁻  for computation of surface processes
        upd_a_surf, upd_e_surf, upd_q_tot_surf =
            compute_updraft_surface_BC(turbconv.surface, turbconv, m, gm, gm_a)
        upd_a_surf = FT(0.1)
        upd_θ_liq_surf = FT(300)
        upd_q_tot_surf = FT(0.0016)
        for i in 1:N
            up[i].ρau = SVector(0, 0, 0)
            up[i].ρa = upd_a_surf
            up[i].ρaθ_liq = upd_θ_liq_surf
            up[i].ρaq_tot = upd_q_tot_surf
        end

        tke, θ_liq_cv, q_tot_cv, θ_liq_q_tot_cv =
            env_surface_covariances(turbconv.surface, turbconv, m, gm, gm_a)
        en_area = environment_area(gm, gm_a, N)
        tke = FT(0)
        θ_liq_cv = FT(0)
        q_tot_cv = FT(0)
        θ_liq_q_tot_cv = FT(0)

        en.ρatke            = gm.ρ * en_area * tke
        en.ρaθ_liq_cv       = gm.ρ * en_area * θ_liq_cv
        en.ρaq_tot_cv       = gm.ρ * en_area * q_tot_cv
        en.ρaθ_liq_q_tot_cv = gm.ρ * en_area * θ_liq_q_tot_cv

    elseif bctype == 2 # top
        ρinv = 1 / gm.ρ
        for i in 1:N
            up[i].ρau = SVector(0, 0, 0)
            up[i].ρa = FT(0)
            up[i].ρaθ_liq = FT(0)
            up[i].ρaq_tot = FT(0)
        end
        en.ρatke = FT(0)
        en.ρaθ_liq_cv = FT(0)
        en.ρaq_tot_cv = FT(0)
        en.ρaθ_liq_q_tot_cv = FT(0)
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
    kernel_calls[:turbconv_normal_boundary_flux_second_order!] = true
    debug_kernels && println("Calling turbconv_normal_boundary_flux_second_order!")
    debug_kernels && @show kernel_calls

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
