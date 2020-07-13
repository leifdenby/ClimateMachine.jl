"""
    DiagnosticVar

A diagnostic variable is an `n`-dimensional field computed from
the various `ClimateMachine` state variables. Diagnostic variables
may be assembled into `DiagnosticsGroup`s.

Standard names are from:
http://cfconventions.org/Data/cf-standard-names/71/build/cf-standard-name-table.html
"""

abstract type DiagnosticVarKind end
struct HorizontalAverage <: DiagnosticVarKind end
struct Diagnostic3D <: DiagnosticVarKind end
struct Diagnostic2D <: DiagnosticVarKind end
struct Variance <: DiagnosticVarKind end
struct Covariance <: DiagnosticVarKind end

struct DiagnosticVar{F}
    kind::DiagnosticVarKind
    name::String
    attrib::OrderedDict
    impl::F

    DiagnosticVar(
        impl,
        kind::DiagnosticVarKind,
        name::String,
        attrib::OrderedDict = OrderedDict(),
    ) = new(kind, name, attrib, impl)
end
const AllDiagnosticVars = OrderedDict{Symbol, DiagnosticVar}()

macro horizontal_average(name, attrib, impl)
    quote
        AllDiagnosticVars[name] =
            DiagnosticVar(HorizontalAverage(), name, $attrib, $impl)
    end
end

@horizontal_average(
    u,
    var_attrib("", "", ""),
    () -> conservative.ρu[1],
)
@horizontal_average(
    w_ht_sgs,
    var_attrib("", "", ""),
    (atmos, conservative, gradient_flux, auxiliary, curr_time) -> begin
        ν, D_t, _ = turbulence_tensors(
            atmos,
            conservative,
            gradient_flux,
            auxiliary,
            curr_time,
        )
        d_h_tot = -D_t .* gradient_flux.∇h_tot
        d_h_tot[end]
    end,
)

@intermediate(
    "thermo",
    (atmos, conservative, auxiliary) -> begin
        thermo_state(atmos, conservative, auxiliary)
    end,
)
@diagnostic3d(
    temp,
    var_attrib("", "", ""),
    (atmos, conservative, auxiliary) -> begin
    end,
)

function vars_atmos_les_default_simple(m::AtmosModel, FT)
    @vars begin
        u::FT
        v::FT
        w::FT
        avg_rho::FT             # ρ
        rho::FT                 # ρρ
        temp::FT
        pres::FT
        thd::FT                 # θ_dry
        et::FT                  # e_tot
        ei::FT                  # e_int
        ht::FT
        hi::FT
        w_ht_sgs::FT

        moisture::vars_atmos_les_default_simple(m.moisture, FT)
    end
end
vars_atmos_les_default_simple(::MoistureModel, FT) = @vars()
function vars_atmos_les_default_simple(m::EquilMoist, FT)
    @vars begin
        qt::FT                  # q_tot
        ql::FT                  # q_liq
        qv::FT                  # q_vap
        thv::FT                 # θ_vir
        thl::FT                 # θ_liq
        w_qt_sgs::FT
    end
end

@diagnostic_group "AtmosLESDefault" vars_atmos_les_default_simple(m, FT) 

