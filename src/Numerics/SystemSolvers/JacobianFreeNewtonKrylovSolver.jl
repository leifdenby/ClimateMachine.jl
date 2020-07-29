
export BatchedJacobianFreeNewtonKrylovSolver

mutable struct JacobianAction{F, FT}
    f!::F
    ϵ::FT
    cache_Fqdq
    cache_Fq
end

"""
Approximations the action of the Jacobian of a nonlinear
form on a vector `Δq` using the difference quotient:

∂F(q)      F(q + ϵΔq) - F(q)
---- Δq ≈ -------------------
 ∂q                ϵ

"""

function (op::JacobianAction)(dQ, Q, args...)
    f! = op.f!
    Fq = op.cache_Fq
    Fqdq = op.cache_Fqdq
    ϵ = op.ϵ

    f!(Fq, Q, args..., increment = false)
    f!(Fqdq, Q + ϵ * dQ, args..., increment = false)
    dQ .= (Fqdq .- Fq) ./ ϵ
end

mutable struct BatchedJacobianFreeNewtonKrylovSolver{ET, AT} <: AbstractNonlinearSolver
    # Tolerance
    ϵ::ET
    # Max number of Newton iterations
    M::Int
    # Maximum number of batched Newton methods (number of vertical columns)
    batch_size::Int
    jvp!
    # Linear solver for the Jacobian system
    linearsolver
    residual::AT
    function BatchedJacobianFreeNewtonKrylovSolver(
        Q::AT,
        jvp!,
        linearsolver,
        batch_size::Int;
        ϵ = √eps(eltype(AT)),
        M = min(20, eltype(Q)),
    ) where {AT}
        residual = similar(Q)
        ET = typeof(ϵ)
        return new{ET, AT}(ϵ, M, batch_size, jvp!, linearsolver, residual)
    end
end

function initialize!(
    implicitoperator!,
    Q,
    Qrhs,
    solver::BatchedJacobianFreeNewtonKrylovSolver,
    args...,
)
    # where R = Qrhs - F(Q)
    R = solver.residual
    # Computes F(Q) and stores in R
    implicitoperator!(R, Q, args...)
    # Computes R = R - Qrhs
    R .-= Qrhs
    return norm(R)
end

function apply_jacobian!(
    JΔQ,
    implicitoperator!,
    Q,
    dQ,
    ϵ,
    args...,
)
    Fq = similar(Q)
    Fqdq = similar(Q)
    implicitoperator!(Fq, Q, args..., increment = false)
    implicitoperator!(Fqdq, Q .+ ϵ .* dQ, args..., increment = false)
    JΔQ .= (Fqdq .- Fq) ./ ϵ
end

function donewtoniteration!(
    implicitoperator!,
    Q,
    Qrhs,
    solver::BatchedJacobianFreeNewtonKrylovSolver,
    tol,
    args...,
)
    ΔQ = similar(Q)
    ϵ = solver.ϵ

    jvp! = ΔQ -> apply_jacobian(implicitoperator!,
        Q,
        ΔQ,
        ϵ,
        args...,
    )

    # Compute right-hand side for Jacobian system:
    # JΔQ = -R
    # where R = Qrhs - F(Q)
    R = solver.residual
    # Computes F(Q) and stores in R
    implicitoperator!(R, Q, args...)
    # Computes R = R - Qrhs
    R .-= Qrhs

    """
    Want:
        (*) linearoperator!(Result, CurrentState:ΔQ, args...)
    
    Want the Jacobian action (jvp!) to behave just like
    a standard rhs evaluation as in (*)
    """

    # Solve JΔQ = -R, here JΔq = F(Q+ϵΔq) - F(Q)/ϵ
    jvp! = (JΔQ, ΔQ, args...) -> apply_jacobian!(
        JΔQ,
        implicitoperator!,
        Q,
        ΔQ,
        ϵ,
        args...,
    )

    linearsolve!(
        jvp!,
        solver.linearsolver,
        ΔQ,
        R,
        args...,
    )

    # Newton correction
    Q .+= ΔQ

    # Reevaluate residual with new solution
    implicitoperator!(R, Q, args...)
    resnorm = norm(R, weighted_norm)
    return resnorm
end
