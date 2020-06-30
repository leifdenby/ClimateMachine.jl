## Courtesy of Charlie:
## https://github.com/CliMA/ClimateMachine.jl/blob/agl/ck/soil_water/tutorials/Land/helper_funcs.jl#L154-L190

using Interpolations

"""
    TimeContinuousData{FT<:AbstractFloat,A<:AbstractArray}
Creates a time-continuous representation
of temporally discrete data. Example:
```julia
FT = Float64
data_discrete = FT[1, 2, 3, 4]
time_discrete = FT[0, 10, 20, 30]
data_continuous = TimeContinuousData(time_discrete, data_discrete)
data_at_40 = data_continuous(40)
```
"""
struct TimeContinuousData{FT<:AbstractFloat,A<:AbstractArray}
  itp
  ext
  bounds::Tuple{FT,FT}
  function TimeContinuousData(time_data::A, data::A) where {A<:AbstractArray}
    FT = eltype(A)
    itp = interpolate((time_data,), data, Gridded(Linear()))
    ext = extrapolate(itp, Flat())
    bounds = (first(data), last(data))
    return new{FT,A}(itp, ext, bounds)
  end
end
function (cont_data::TimeContinuousData)(x::A) where {A<:AbstractArray}
  return [ cont_data.bounds[1] < x_i < cont_data.bounds[2] ?
  cont_data.itp(x_i) : cont_data.ext(x_i) for x_i in x]
end

function (cont_data::TimeContinuousData)(x_i::FT) where {FT<:AbstractFloat}
  return cont_data.bounds[1] < x_i < cont_data.bounds[2] ? cont_data.itp(x_i) : cont_data.ext(x_i)
end
