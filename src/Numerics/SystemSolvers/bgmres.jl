#### Generalized Minimal Residual Solver

export BatchedGeneralizedMinimalResidual

mutable struct BatchedGeneralizedMinimalResidual{M, MP1, MMP1, T, AT} <:
               AbstractIterativeSystemSolver
    krylov_basis::NTuple{MP1, AT}
    "Hessenberg matrix"
    H::Vector{Matrix{T}}
    "rhs of the least squares problem"
    g0::Vector{Vector{T}}
    rtol::T
    atol::T
    batch_size::Int
    forward_reshape
    forward_permute

    function BatchedGeneralizedMinimalResidual(
        Q::AT,
        Nbatch;
        M = min(20, eltype(Q)),
        rtol = √eps(eltype(AT)),
        atol = eps(eltype(AT)),
        forward_reshape = size(Q),
        forward_permute = Tuple(1:length(size(Q))),
    ) where {AT}

        #(k_n + 1,  m,  n)
        krylov_basis = ntuple(i -> similar(Q), M + 1)
        H = zeros(Nbatch, M + 1, M)
        g0 = zeros(Nbatch, M + 1)

        new{M, M + 1, M * (M + 1), eltype(Q), AT}(
            krylov_basis,
            H,
            g0,
            rtol,
            atol,
            Nbatch,
            forward_reshape,
            forward_permute,
        )
    end
end

"""
    BatchedGeneralizedMinimalResidual(
        dg::DGModel,
        Q::MPIStateArray;
        atol = sqrt(eps(eltype(Q))),
        rtol = sqrt(eps(eltype(Q))),
        max_subspace_size = nothing,
        independent_states = false,
    )

# Description
Specialized constructor for `BatchedGeneralizedMinimalResidual` struct, using
a `DGModel` to infer state-information and determine appropriate reshaping
and permutations.

# Arguments
- `dg`: (DGModel) A `DGModel` containing all relevant grid and topology
        information.
- `Q` : (MPIStateArray) An `MPIStateArray` containing field information.

# Keyword Arguments
- `atol`: (float) absolute tolerance. `DEFAULT = sqrt(eps(eltype(Q)))`
- `rtol`: (float) relative tolerance. `DEFAULT = sqrt(eps(eltype(Q)))`
- `max_subspace_size` : (Int).    Maximal dimension of each (batched)
                                  Krylov subspace. DEFAULT = nothing
- `independent_states`: (boolean) An optional flag indicating whether
                                  or not degrees of freedom are coupled
                                  internally (within a column).
                                  `DEFAULT = false`
# Return
instance of `BatchedGeneralizedMinimalResidual` struct
"""
function BatchedGeneralizedMinimalResidual(
    dg::DGModel,
    Q::MPIStateArray;
    atol = sqrt(eps(eltype(Q))),
    rtol = sqrt(eps(eltype(Q))),
    max_subspace_size = nothing,
    independent_states = false,
)

    # Need to determine array type for storage vectors
    if isa(Q.data, Array)
        ArrayType = Array
    else
        ArrayType = CuArray
    end

    grid = dg.grid
    topology = grid.topology
    dim = dimensionality(grid)

    # Number of Gauss-Lobatto quadrature points in 1D
    Nq = polynomialorder(grid) + 1

    # Assumes same number of quadrature points in all spatial directions
    Np = Tuple([Nq for i in 1:dim])

    # Number of states and elements (in vertical and horizontal directions)
    num_states = size(Q)[2]
    nelem = length(topology.elems)
    nvertelem = topology.stacksize
    nhorzelem = div(nelem, nvertelem)

    # Definition of a "column" here is a vertical stack of degrees
    # of freedom. For example, consider a mesh consisting of a single
    # linear element:
    #    o----------o
    #    |\ d1   d2 |\
    #    | \        | \
    #    |  \ d3    d4 \
    #    |   o----------o
    #    o--d5---d6-o   |
    #     \  |       \  |
    #      \ |        \ |
    #       \|d7    d8 \|
    #        o----------o
    # There are 4 total 1-D columns, each containing two
    # degrees of freedom. In general, a mesh of stacked elements will
    # have `Nq^2 * nhorzelem` total 1-D columns.
    # A single 1-D column has `Nq * nvertelem * num_states`
    # degrees of freedom.
    #
    # nql = length(Np)
    # indices:      (1...nql, nql + 1 , nql + 2, nql + 3)

    # for 3d case, this is [ni, nj, nk, num_states, nvertelem, nhorzelem]
    # here ni, nj, nk are number of Gauss quadrature points in each element in x-y-z directions
    # Q = reshape(Q, reshaping_tup), leads to the column-wise fashion Q
    reshaping_tup = (Np..., num_states, nvertelem, nhorzelem)

    if independent_states
        m = Nq * nvertelem
        n = (Nq^(dim - 1)) * nhorzelem * num_states
    else
        m = Nq * nvertelem * num_states
        n = (Nq^(dim - 1)) * nhorzelem
    end

    if max_subspace_size === nothing
        max_subspace_size = m
    end


    # permute [ni, nj, nk, num_states, nvertelem, nhorzelem] 
    # to      [nvertelem, nk, num_states, ni, nj, nhorzelem]
    permute_size = length(reshaping_tup)
    permute_tuple_f = (dim + 2, dim, dim + 1, (1:dim-1)..., permute_size)

    # Now we need to determine an appropriate permutation
    # of the MPIStateArray to perform column-wise strides.
    # Total size of the permute tuple
    permute_size = length(reshaping_tup)

    # Index associated with number of GL points
    # in the 'vertical' direction
    nql = length(Np)
    # Want: (index associated with the stack size of the column,
    #        index associated with GL pts in vertical direction,
    #        index associated with number of states)
    # FIXME: Better way to do this?
    #             (vert stack, GL pts, num states)
    column_strides = (nql + 2, nql, nql + 1)
    diff = Tuple(setdiff(Set([i for i in 1:permute_size]), Set(column_strides)))
    permute_tuple_f = (column_strides..., diff...)

    return BatchedGeneralizedMinimalResidual(
        Q;
        m = m,
        n = n,
        subspace_size = max_subspace_size,
        atol = atol,
        rtol = rtol,
        ArrayType = ArrayType,
        reshape_tuple_f = reshaping_tup,
        permute_tuple_f = permute_tuple_f,
    )
end



function initialize!(
    linearoperator!,
    Q,
    Qrhs,
    solver::GeneralizedMinimalResidual,
    args...,
)
    g0 = solver.g0
    krylov_basis = solver.krylov_basis
    rtol, atol = solver.rtol, solver.atol

    @assert size(Q) == size(krylov_basis[1])

    # store the initial residual in krylov_basis[1]
    linearoperator!(krylov_basis[1], Q, args...)
    @. krylov_basis[1] = Qrhs - krylov_basis[1]

    threshold = rtol * norm(krylov_basis[1], weighted_norm)
    residual_norm = norm(krylov_basis[1], weighted_norm)

    converged = false
    # FIXME: Should only be true for threshold zero
    if threshold < atol
        converged = true
        return converged, threshold
    end

    fill!(g0, 0)
    g0[1] = residual_norm
    krylov_basis[1] ./= residual_norm

    converged, max(threshold, atol)
end

function doiteration!(
    linearoperator!,
    Q,
    Qrhs,
    solver::GeneralizedMinimalResidual{M},
    threshold,
    args...,
) where {M}

    krylov_basis = solver.krylov_basis
    H = solver.H
    g0 = solver.g0

    converged = false
    residual_norm = typemax(eltype(Q))

    Ω = LinearAlgebra.Rotation{eltype(Q)}([])
    j = 1
    for outer j in 1:M

        # Arnoldi using the Modified Gram Schmidt orthonormalization
        linearoperator!(krylov_basis[j + 1], krylov_basis[j], args...)
        for i in 1:j
            H[i, j] = dot(krylov_basis[j + 1], krylov_basis[i], weighted_norm)
            @. krylov_basis[j + 1] -= H[i, j] * krylov_basis[i]
        end
        H[j + 1, j] = norm(krylov_basis[j + 1], weighted_norm)
        krylov_basis[j + 1] ./= H[j + 1, j]

        # apply the previous Givens rotations to the new column of H
        @views H[1:j, j:j] .= Ω * H[1:j, j:j]

        # compute a new Givens rotation to zero out H[j + 1, j]
        G, _ = givens(H, j, j + 1, j)

        # apply the new rotation to H and the rhs
        H .= G * H
        g0 .= G * g0

        # compose the new rotation with the others
        Ω = lmul!(G, Ω)

        residual_norm = abs(g0[j + 1])

        if residual_norm < threshold
            converged = true
            break
        end
    end

    # solve the triangular system
    y = SVector{j}(@views UpperTriangular(H[1:j, 1:j]) \ g0[1:j])

    ## compose the solution
    rv_Q = realview(Q)
    rv_krylov_basis = realview.(krylov_basis)
    groupsize = 256
    event = Event(array_device(Q))
    event = linearcombination!(array_device(Q), groupsize)(
        rv_Q,
        y,
        rv_krylov_basis,
        true;
        ndrange = length(rv_Q),
        dependencies = (event,),
    )
    wait(array_device(Q), event)

    # if not converged restart
    converged || initialize!(linearoperator!, Q, Qrhs, solver, args...)

    (converged, j, residual_norm)
end
