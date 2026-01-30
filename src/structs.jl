
struct SplitTime{T<:Union{Int, Vector{<:Real}}}
    spec::T
end
"""
    SplitTime(n::Integer)

Split the measurements in the measurement table into `n` chunks based on unique time points.

The resulting chunking is applied to each simulation condition. If the number of unique time
points is not divisible by `n`, the last chunk contains the remaining time points.
"""
SplitTime(n::T) where T <: Integer = SplitTime{T}(n)
"""
    SplitTime(split_points::Vector{<:Real})

Splits the measurements in the measurement table into `length(split_points) + 1` chunks
using the time values in `split_points` as split points.

Chunks are defined as  `[tmin, s₁)`, `[s₁, s₂)`, …, `[sₙ, tmax]`, where `tmin`/`tmax` are
the minimum/maximum measurement time points and `sᵢ` are the entries in `split_points`.
`split_points` must be sorted and lie within the range of the measurement time points.
Split points do not need to coincide with measurement time points.
"""
SplitTime(split_points::T) where T <: Vector{<:Real} = SplitTime{T}(split_points)

struct SplitData{T<:Union{Int, Vector{<:Integer}}}
    spec::T
end
"""
    SplitData(n::Integer)

Splits the time-sorted measurements into `n` chunks (uniform by number of data points).

The resulting chunking is applied to each simulation condition. If the number of data points
is not divisible by `n`, the last chunk contains the remaining data points.

If all measurements have a unique time point, this is equivalent to `SplitTime(n)`.
"""
SplitData(n::T) where T <: Integer = SplitData{T}(n)
"""
    SplitData(chunk_sizes::Vector{<:Integer})

Splits the measurements in the time-sorted measurement table into `length(chunk_sizes)`
chunks, where `chunk_sizes[i]` gives the number of data points in chunk `i`.

`chunk_sizes` must be sorted, and the final entry must equal the number of measurements.
"""
SplitData(chunk_sizes::T) where T <: Vector{<:Integer} = SplitData{T}(chunk_sizes)

"""
    MsInitConstant(value::Real = 0.01)

Initialization strategy for multiple-shooting window initial values, where each window
initial value is set to `value`.

See also [`set_u0_ms_windows!`](@ref).
"""
struct MsInitConstant{T<:Real}
    value::T
end
MsInitConstant() = MsInitConstant{Float64}(0.01)

"""
    MsInitFirst()

Initialization strategy for multiple-shooting window initial values, where all windows
take the initial values of the first window (values at `t0`).

See also [`set_u0_ms_windows!`](@ref).
"""
struct MsInitFirst end

"""
    MsInitSimulate()

Initialization strategy for multiple-shooting window initial values, where each window
initial value is set from a forward simulation of the model up to the start time of the
corresponding window.

The model is simulated using options stored in the PEtab problem.

See also [`set_u0_ms_windows!`](@ref).
"""
struct MsInitSimulate end
