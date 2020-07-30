
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
    # nonlinear rhs
    rhs!
    # Linear solver for the Jacobian system
    linearsolver
    residual::AT
    function BatchedJacobianFreeNewtonKrylovSolver(
        Q::AT,
        rhs!,
        linearsolver,
        batch_size::Int;
        ϵ = √eps(eltype(AT)),
        M = min(20, eltype(Q)),
    ) where {AT}
        residual = similar(Q)
        ET = typeof(ϵ)
        return new{ET, AT}(ϵ, M, batch_size, rhs!, linearsolver, residual)
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

function donewtoniteration!(
    implicitoperator!,
    jvp!,
    Q,
    Qrhs,
    solver::BatchedJacobianFreeNewtonKrylovSolver,
    tol,
    args...,
)
    FT = eltype(Q)
    ΔQ = similar(Q)
    ΔQ .= FT(0.0)

    # Compute right-hand side for Jacobian system:
    # JΔQ = -R
    # where R = Qrhs - F(Q)
    R = solver.residual
    # Computes F(Q) and stores in R
    implicitoperator!(R, Q, args...)
    # Computes R = R - Qrhs
    R .-= Qrhs

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
