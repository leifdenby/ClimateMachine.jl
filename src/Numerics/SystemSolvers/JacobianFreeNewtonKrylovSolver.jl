
export BatchedJacobianFreeNewtonKrylovSolver

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
    return false, zero(eltype(Q))
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
    jvp! = solver.jvp!

    # Compute right-hand side for Jacobian system:
    # JΔQ = -R
    # where R = Qrhs - F(Q)
    R = solver.residual
    implicitoperator!(R, Q, args...)
    R .-= Qrhs

    # Solve JΔq = -R
    linearsolve!(
        jvp!,
        solver.linearsolver,
        ΔQ,
        R,
        args...;
        max_iters = getmaxiterations(solver.linearsolver),
    )

    # Newton correction
    Q .+= ΔQ

    # Reevaluate residual with new solution
    implicitoperator!(R, Q, args...)
    resnorm = norm(R, weighted_norm)
    converged = false
    if resnorm < tol
        converged = true
    end

    return converged, resnorm
end
