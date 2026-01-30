struct SplitTime{T <: Union{Int, Vector{<:Real}}}
    spec::T
end
"""
    SplitTime(n::Integer)

Split the measurements in the measurement table into `n` chunks based on unique time points.

The resulting chunking is applied to each simulation condition. If the number of unique time
points is not divisible by `n`, the last chunk contains the remaining time points.
"""
SplitTime(n::T) where {T <: Integer} = SplitTime{T}(n)
"""
    SplitTime(split_points::Vector{<:Real})

Splits the measurements in the measurement table into `length(split_points) + 1` chunks
using the time values in `split_points` as split points.

Chunks are defined as  `[tmin, s₁)`, `[s₁, s₂)`, …, `[sₙ, tmax]`, where `tmin`/`tmax` are
the minimum/maximum measurement time points and `sᵢ` are the entries in `split_points`.
`split_points` must be sorted and lie within the range of the measurement time points.
Split points do not need to coincide with measurement time points.
"""
SplitTime(split_points::T) where {T <: Vector{<:Real}} = SplitTime{T}(split_points)

struct SplitData{T <: Union{Int, Vector{<:Integer}}}
    spec::T
end
"""
    SplitData(n::Integer)

Splits the time-sorted measurements into `n` chunks (uniform by number of data points).

The resulting chunking is applied to each simulation condition. If the number of data points
is not divisible by `n`, the last chunk contains the remaining data points.

If all measurements have a unique time point, this is equivalent to `SplitTime(n)`.
"""
SplitData(n::T) where {T <: Integer} = SplitData{T}(n)
"""
    SplitData(chunk_sizes::Vector{<:Integer})

Splits the measurements in the time-sorted measurement table into `length(chunk_sizes)`
chunks, where `chunk_sizes[i]` gives the number of data points in chunk `i`.

`chunk_sizes` must be sorted, and the final entry must equal the number of measurements.
"""
SplitData(chunk_sizes::T) where {T <: Vector{<:Integer}} = SplitData{T}(chunk_sizes)

"""
    MsInitConstant(value::Real = 0.01)

Initialization strategy for multiple-shooting window initial values, where each window
initial value is set to `value`.

See also [`set_u0_ms_windows!`](@ref).
"""
struct MsInitConstant{T <: Real}
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

"""
    PEtabClProblem(prob_original::PEtabODEProblem, split_alg::Union{SplitTime, SplitData})

Construct a curriculum of PEtab problems by splitting the measurement data of
`prob_original` into curriculum stages using `split_alg`.

Stage problems are stored in `petab_problems`, where `petab_problems[i]` is the
`PEtabODEProblem` for curriculum stage `i`. The original problem is stored in `original`.

See also [`SplitTime`](@ref), [`SplitData`](@ref).
"""
struct PEtabClProblem{T <: Union{SplitTime, SplitData}}
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
    split_alg::T
    regularization_obs::Symbol
end

"""
    PEtabMsProblem(prob_original::PEtabODEProblem, split_alg::SplitTime; kwargs...)

Construct a multiple-shooting PEtab problem by splitting the time domain of `prob_original`
into windows defined by `split_alg`.

In the constructed `prob_ms`, the multiple-shooting `PEtabODEProblem` is stored in
`prob_ms.petab_ms_problem`, and the original problem is stored in `prob_ms.original`.

# Keyword arguments
- `window_u0_scale::Symbol = :lin`: Scale used for the estimated window initial-value
  parameters (`:lin`, `:log`, `:log10`). For problems where ODE states span orders
  of  magnitude, `:log10` often improves performance. Any other option than `:lin`
  should only be applied if states are positive as it enforces positive state-values.

# Mathematical description

Multiple shooting introduces additional window initial-value parameters that are estimated
together with the original parameters. Window continuity is enforced by adding a quadratic
penalty at each window boundary. For consecutive windows with end time `t_i` and next-window
initial-value parameter `\tilde{u}_{i+1}`, the penalty used is

```math
\\lambda \\left\\lVert u_i(t_i) - \\tilde{u}_{i+1} \\right\\rVert_2^2,
```

where ``u_i(t_i)`` is the state at the end of window ``i``, and ``\\lambda`` is the
window-penalty weight (can be set with [`set_window_penalty!`](@ref)).

Initial values for the window parameters `\tilde{u}_{i+1}` can be set with
[`set_u0_ms_windows!`](@ref) (default `0.01`).

See also [`SplitTime`](@ref), [`set_window_penalty!`](@ref), [`set_u0_ms_windows!`](@ref).
"""
struct PEtabMsProblem
    petab_ms_problem::PEtabODEProblem
    original::PEtabODEProblem
    split_alg::SplitTime
    ms_windows::Vector{Vector{Float64}}
    window_u0_scale::Symbol
end

"""
    PEtabClMsProblem(prob_original::PEtabODEProblem, split_alg::SplitTime; kwargs...)

Construct a curriculum of multiple-shooting PEtab problems from `prob_original` by
progressively merging time windows defined by `split_alg`.

In the constructed `prob_cl_ms`, stage problems are stored in `prob_cl_ms.petab_problems`,
where `prob_cl_ms.petab_problems[i]` is the `PEtabODEProblem` for curriculum stage `i`. The
original problem is stored in `prob_cl_ms.original`.

Each intermediate stage is a multiple-shooting problem with a different window partition
(with increasing window length and overlap). The final stage equals the original problem.

As in [`PEtabMsProblem`](@ref), each multiple-shooting stage includes a quadratic window
penalty applied at the first time point shared by two consecutive windows; the penalty
weight can be set with [`set_window_penalty!`](@ref).

Because the number of shooting windows changes between stages, the parameter vector dimension
changes as well. Use [`map_x_stage`](@ref) to map a parameter vector between stages.

Keyword arguments are forwarded to [`PEtabMsProblem`](@ref) (`window_u0_scale`).

See also [`SplitTime`](@ref), [`PEtabMsProblem`](@ref), [`set_window_penalty!`](@ref),
[`map_x_stage`](@ref).
"""
struct PEtabClMsProblem
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
    split_alg::SplitTime
    ms_windows::Dict{Symbol, Vector{Vector{Float64}}}
    window_u0_scale::Symbol
end
