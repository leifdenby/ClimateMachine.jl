dofs_per_elem_per_reduced_direction(N::Int, dim::Int) =
    dim == 2 ? 1 : dofs_per_elem_per_direction(N) # better name for this method?

dofs_per_elem_reduced(N::Int, dim::Int) = # better name for this method?
    dofs_per_elem_per_direction(N) * dofs_per_elem_per_reduced_direction(N, dim)

dofs_per_elem_per_direction(N::Int) = (N + 1) # better name for this method?

dofs_per_elem(N::Int, dim::Int) = (N + 1)^dim

n_faces(dim::Int) = 2 * dim

get_faces(dim::Int, direction) = 1:n_faces(dim)

get_faces(dim::Int, ::VerticalDirection) = (n_faces(dim) - 1):n_faces(dim)

get_faces(dim::Int, ::HorizontalDirection) = 1:(n_faces(dim) - 2)

"""
    get_face_direction(face, dim)

The remainder model needs to know which direction of face the model is
being evaluated for. So faces 1:(n_faces(dim) - 2) are flagged as
`HorizontalDirection()` faces and the remaining two faces are
`VerticalDirection()` faces
"""
function get_face_direction(face::Int, dim::Int)
    if face in 1:(n_faces(dim) - 2)
        return (EveryDirection(), HorizontalDirection())
    else
        return (EveryDirection(), VerticalDirection())
    end
end
