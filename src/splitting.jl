struct SplitUniform
    nsplits::Integer
    mode::Symbol
end
function SplitUniform(nsplits::Integer; mode::Symbol = :time)::SplitUniform
    @argcheck mode in [:time, :condition] "Mode must be either :time or :condition"
    return SplitUniform(nsplits, mode)
end
struct SplitCustom{T <: Union{AbstractFloat, String, Symbol}}
    nsplits::Integer
    mode::Symbol
    splits::Vector{Vector{T}}
end
function SplitCustom(splits::Union{Vector, Vector{<:Vector}})::SplitCustom
    # For curriculum the only information that must be provided is the end time-point
    # for each curriculum step. Internally, for a common interface, splits must be a
    # Vector{Vector}.
    if splits isa Vector
        splits = [[0.0, splits[i]] for i in eachindex(splits)]
    end
    if !(all(eltype.(splits) .<: Real) || all(eltype.(splits) .<: String) || all(eltype.(splits) .<: Symbol))
        throw(ArgumentError("Provided splits must be a Vector{Vector} where the inner \
            Vector all are of the same type; either Vector{Real} if splitting on time, \
            or either Vector{String} or Vector{Symbol} if splitting on conditions"))
    end
    for (i, split) in pairs(splits)
        !isempty(split) && continue
        throw(ArgumentError("Provided interval $i is empty. When providing custom \
            splitting intervals no interval must be empty"))
    end
    nsplits = length(splits)
    mode = all(eltype.(splits) .<: Real) ? :time : :condition
    return SplitCustom(nsplits, mode, splits)
end

function _split(split_algorithm::SplitUniform, prob::PEtabODEProblem)
    if split_algorithm.mode == :time
        return _split_uniform_time(prob, split_algorithm.nsplits)
    end
    return _split_uniform_conditions(prob, split_algorithm.nsplits)
end
function _split(split_algorithm::SplitCustom, prob::PEtabODEProblem)
    if split_algorithm.mode == :time
        return _split_custom_time(prob, split_algorithm.splits, split_algorithm.nsplits)
    end
    return _split_custom_conditions(prob, split_algorithm.splits, split_algorithm.nsplits)
end

function _split_uniform_time(prob::PEtabODEProblem, nsplits::Integer)
    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_timepoints(mdf)
    if length(unique_t) < nsplits
        throw(ArgumentError("Number of unique time-points in the measurement table must be \
          greater than or equal to number of splits (number of curriculum stages or \
          multiple shooting windows). There are $(length(unique_t)) unique \
          time-points and nsplits=$(nsplits)."),
        )
    end
    splits = _makechunks(unique_t, nsplits)
    return _split_time(splits, mdf, prob)
end

function _split_custom_time(prob::PEtabODEProblem, splits::Vector{<:Vector{<:Real}}, nsplits::Integer)
    for i in 1:(nsplits - 1)
        splits[i][end] < splits[i+1][end] && continue
        throw(ArgumentError("When providing custom time-points for curriculum learning, \
            end point for interval i+1 must be larger than for interval i. This does \
            not hold for i=$i; its end-point is $(splits[i][end]) and the end point for \
            i=$(i+1) is $(splits[i+1][end])"))
    end
    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_timepoints(mdf)
    n_time_points_included = 0
    for (i, split) in pairs(splits)
        tmax_split= maximum(split)
        n_time_points_split = sum(unique_t .≤ tmax_split)
        if n_time_points_split > n_time_points_included
            n_time_points_included = n_time_points_split
            continue
        end
        if i == 1
            throw(ArgumentError("When providing custom time-points for curriculum \
                learning, each curriculum step i must include new measurements points. \
                This does not hold for step $i. This means there are no measurement points \
                in the measurement table between [0.0, $(tmax_split)]."))
        end
        tmax_prev = maximum(splits[i-1])
        throw(ArgumentError("When providing custom time-points for curriculum \
            learning, each curriculum step i must include new measurements points. \
            This does not hold for step $i. This means there are no measurement points \
            in the measurement table between ($(tmax_prev), $(tmax_split)]."))
    end
    if maximum(splits[end]) < maximum(unique_t)
        throw(ArgumentError("When providing custom time-points for curriculum, the last \
            stage must cover all measurements in the measurement table. This does not hold \
            as the end-point time for the last interval is $(maximum(splits[end])) while \
            the largest time-point in the measurements table is $(maximum(unique_t))."))
    end
    return _split_time(splits, mdf, prob)
end

function _split_uniform_conditions(prob::PEtabODEProblem, nsplits::Integer)
    conditions_ids = prob.model_info.simulation_info.conditionids[:simulation] .|> string
    if length(conditions_ids) < nsplits
        throw(ArgumentError("Number of conditions in the measurement table must be \
            greater than or equal to number of splits (number of curriculum stages or \
            multiple shooting windows). There model has $(length(conditions_ids)) unique \
            conditions Ids and nsplits=$(nsplits)."))
    end

    mdf = prob.model_info.model.petab_tables[:measurements]
    conditions_stage = _makechunks(conditions_ids, nsplits)
    out = Vector{Dict{Symbol, DataFrame}}(undef, nsplits)
    for i in 1:nsplits
        out[i] = copy(prob.model_info.model.petab_tables)
        irow = findall(x -> x in conditions_stage[i], mdf.simulationConditionId)
        out[i][:measurements] = mdf[irow, :]
    end
    return out
end

function _get_unique_timepoints(mdf::DataFrame)::Vector{Float64}
    return mdf.time |>
        unique |>
        sort
end

function _split_time(splits, mdf::DataFrame, prob::PEtabODEProblem)::Vector{Dict{Symbol, DataFrame}}
    out = Vector{Dict{Symbol, DataFrame}}(undef, length(splits))
    for (i, split) in pairs(splits)
        tmax = maximum(split)
        out[i] = copy(prob.model_info.model.petab_tables)
        out[i][:measurements] = mdf[mdf[!, :time] .≤ tmax, :]
    end
    return out
end
