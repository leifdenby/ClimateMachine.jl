function environment_area(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    return 1 - sum([state.turbconv.updraft[i].ρa for i in 1:N])/ state.ρ
end

function environment_w(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    return (state.ρu[3] - sum([state.turbconv.updraft[i].ρaw for i in 1:N]))/a_en*ρinv
end

function grid_mean_b(
    state::Vars,
    aux::Vars,
    N::Int,
) where {FT}
    ρinv = 1/state.ρ
    a_en = environment_area(state ,aux ,N)
    up = state.turbconv.updraft
    en_a = aux.turbconv.environment
    up_a = aux.turbconv.updraft
    return a_en * en_a.buoyancy + sum([up_a[i].buoyancy*up[i].ρa*ρinv for i in 1:N])
end