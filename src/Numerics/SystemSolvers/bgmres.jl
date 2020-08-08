
export BatchedGeneralizedMinimalResidual

mutable struct BatchedGeneralizedMinimalResidual{
    MP1,
    I,
    T,
    AT,
    FRS,
    FPR,
    BRS,
    BPR
} <: AbstractIterativeSystemSolver

    "global Krylov basis"
    krylov_basis::NTuple{MP1, AT}
    "Hessenberg matrix in each column"
    H::Vector{Matrix{T}}
    "rhs of the least squares problem in each column"
    g0::Vector{Vector{T}}
    "global work vector for computing A*v"
    globwvec::AT
    rtol::T
    atol::T
    "total number of batched columns"
    batch_size::I
    "iterations to reach convergence in each column"
    iterconv::Vector{I}
    "residual norm in each column"
    resnorms::Vector{T}
    forward_reshape::FRS
    forward_permute::FPR
    backward_reshape::BRS
    backward_permute::BPR

    function BatchedGeneralizedMinimalResidual(
        Q::AT,
        Nbatch,
        dofperbatch;
        M = min(20, length(Q)),
        rtol = √eps(eltype(AT)),
        atol = eps(eltype(AT)),
        forward_reshape = size(Q),
        forward_permute = Tuple(1:length(size(Q))),
    ) where {AT}
        FRS = typeof(forward_reshape)
        FPR = typeof(forward_permute)
        krylov_basis = ntuple(i -> similar(Q), M + 1)
        H = -zeros(eltype(AT), (Nbatch, M + 1, M))
        g0 = -zeros(eltype(AT), (Nbatch, M + 1))

        @assert dofperbatch * Nbatch == length(Q)

        globwvec = -zeros(eltype(AT), (dofperbatch, Nbatch))
        iterconv = fill(-1, Nbatch)
        resnorms = -zeros(eltype(AT), Nbatch)

        # define the back permutation and reshape
        backward_permute = Tuple(sortperm([forward_permute...]))
        tmp_reshape_tuple_b = [forward_reshape...]
        permute!(tmp_reshape_tuple_b, [forward_permute...])
        backward_reshape = Tuple(tmp_reshape_tuple_b)
        BRS = typeof(backward_reshape)
        BPR = typeof(backward_permute)

        new{M + 1, typeof(Nbatch), eltype(Q), AT, FRS, FPR, BRS, BPR}(
            krylov_basis,
            H,
            g0,
            globwvec,
            rtol,
            atol,
            Nbatch,
            iterconv,
            resnorms,
            forward_reshape,
            forward_permute,
            backward_reshape,
            backward_permute,
        )
    end
end

"""
    BatchedGeneralizedMinimalResidual(
        dg::DGModel,
        Q::MPIStateArray;
        atol = sqrt(eps(eltype(Q))),
        rtol = sqrt(eps(eltype(Q))),
        max_iteration = nothing,
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
- `max_iteration` : (Int).    Maximal dimension of each (batched)
                              Krylov subspace. DEFAULT = nothing
# Return
instance of `BatchedGeneralizedMinimalResidual` struct
"""
function BatchedGeneralizedMinimalResidual(
    dg::DGModel,
    Q::MPIStateArray;
    atol = sqrt(eps(eltype(Q))),
    rtol = sqrt(eps(eltype(Q))),
    max_iteration = nothing,
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

    m = Nq * nvertelem * num_states
    n = (Nq^(dim - 1)) * nhorzelem

    if max_iteration === nothing
        max_iteration = m
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
        Q,
        n,
        m;
        M = max_iteration,
        atol = atol,
        rtol = rtol,
        forward_reshape = reshaping_tup,
        forward_permute = permute_tuple_f,
    )
end

@inline function convert_structure!(
    x,
    y,
    reshape_tuple,
    permute_tuple,
)
    alias_y = reshape(y, reshape_tuple)
    permute_y = permutedims(alias_y, permute_tuple)
    x[:] .= permute_y[:]
    nothing
end
@inline convert_structure!(x, y::MPIStateArray, reshape_tuple, permute_tuple) =
    convert_structure!(x, y.data, reshape_tuple, permute_tuple)
@inline convert_structure!(x::MPIStateArray, y, reshape_tuple, permute_tuple) =
    convert_structure!(x.data, y, reshape_tuple, permute_tuple)

function initialize!(
    linearoperator!,
    Q,
    Qrhs,
    solver::BatchedGeneralizedMinimalResidual,
    args...,
)
    g0 = solver.g0
    krylov_basis = solver.krylov_basis
    rtol, atol = solver.rtol, solver.atol

    iterconv = solver.iterconv
    globwvec = solver.globwvec
    forward_reshape = solver.forward_reshape
    forward_permute = solver.forward_permute
    resnorms = solver.resnorms

    # Get device and groupsize information
    device = array_device(Q)
    if isa(device, CPU)
        groupsize = Threads.nthreads()
    else # isa(device, CUDADevice)
        groupsize = 256
    end

    @assert size(Q) == size(krylov_basis[1])

    # FIXME: Can we make linearoperator! batch-able?
    # store the initial residual in krylov_basis[1]
    linearoperator!(krylov_basis[1], Q, args...)
    @. krylov_basis[1] = Qrhs - krylov_basis[1]

    convert_structure!( 
        globwvec,
        krylov_basis[1],
        forward_reshape,
        forward_permute,
    )

    event = Event(device)
    event = batched_initialize!(device, groupsize)(
        resnorms,
        iterconv,
        g0,
        globwvec;
        ndrange = solver.batch_size,
        dependencies = (event,),
    )
    wait(device, event)

    residual_norm = maximum(resnorms)
    threshold = rtol * residual_norm
    converged = false
    if threshold < atol
        converged = true
    end

    converged, max(threshold, atol)
end

@kernel function batched_initialize!(resnorms, iterconv, g0, krylov_basis)
    cidx = @index(Global)

    FT = eltype(krylov_basis[cidx])
    fill!(g0[cidx], FT(0.0))

    local_residual_norm = norm(krylov_basis[cidx], false)
    g0[cidx][1] = local_residual_norm
    @. krylov_basis[cidx] /= local_residual_norm
    resnorms[cidx] = local_residual_norm
    iterconv[cidx] = 0

    nothing
end
