struct SplitUniform
    nsplits::Integer
    mode::Symbol
end
function SplitUniform(nsplits::Integer; mode::Symbol = :time)::SplitUniform
    @argcheck mode in [:time, :condition, :datapoints] "Mode must be either :time, \
         :datapoints or :condition"
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
    if splits isa Vector{<:Real}
        splits = [[0.0, splits[i]] for i in eachindex(splits)]
    end
    if !(all(eltype.(splits) .<: Real) || all(eltype.(splits) .<: String) ||
         all(eltype.(splits) .<: Symbol))
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

function _split(split_algorithm::SplitUniform, prob::PEtabODEProblem, method::Symbol)
    mode = split_algorithm.mode
    if method == :multiple_shooting && (mode == :condition || mode == :datapoints)
        throw(ArgumentError("For multiple shotting, splitting windows over $mode \
            (mode = :$(mode)) is not allowed."))
    end

    if mode == :time
        return _split_uniform_time(prob, split_algorithm.nsplits, method)
    elseif mode == :condition
        return _split_uniform_conditions(prob, split_algorithm.nsplits)
    elseif mode == :datapoints
        return _split_uniform_datapoints(prob, split_algorithm.nsplits)
    end
end
function _split(split_algorithm::SplitCustom, prob::PEtabODEProblem, method::Symbol)
    if split_algorithm.mode == :time
        return _split_custom_time(prob, split_algorithm.splits, split_algorithm.nsplits,
            method::Symbol)
    end
    if method == :multiple_shooting
        throw(ArgumentError("For multiple shotting, splitting windows over conditions \
            (mode = :condition) is not allowed."))
    end
    return _split_custom_conditions(prob, split_algorithm.splits)
end

function _split_uniform_time(prob::PEtabODEProblem, nsplits::Integer, method::Symbol)
    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_timepoints(mdf)
    if length(unique_t) < nsplits
        throw(ArgumentError("Number of unique time-points in the measurement table must be \
          greater than or equal to number of splits (number of curriculum stages or \
          multiple shooting windows). There are $(length(unique_t)) unique \
          time-points and nsplits=$(nsplits)."),
        )
    end
    if method == :multiple_shooting
        return _makechunks(unique_t, nsplits; overlap = 1)
    else
        splits = _makechunks(unique_t, nsplits; overlap = 0)
        return _split_time_curriculum(splits, mdf, prob)
    end
end

function _split_custom_time(prob::PEtabODEProblem, splits::Vector{<:Vector{<:Real}},
        nsplits::Integer, method::Symbol)
    for i in 1:(nsplits - 1)
        splits[i][end] < splits[i + 1][end] && continue
        throw(ArgumentError("When providing custom time-points for curriculum learning or \
            multiple-shooting,end point for interval i+1 must be larger than for interval \
            i. This does not hold for i=$i; its end-point is $(splits[i][end]) and the
            end point for i=$(i+1) is $(splits[i+1][end])"))
    end
    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_timepoints(mdf)
    n_time_points_included = 0
    for (i, split) in pairs(splits)
        tmax_split = maximum(split)
        n_time_points_split = sum(unique_t .≤ tmax_split)
        if n_time_points_split > n_time_points_included
            n_time_points_included = n_time_points_split
            continue
        end
        if i == 1
            throw(ArgumentError("When providing custom time-points for curriculum learning \
                or multiple-shooting, step/window i must include new measurements points. \
                This does not hold for step $i. This means there are no measurement points \
                in the measurement table between [0.0, $(tmax_split)]."))
        end
        tmax_prev = maximum(splits[i - 1])
        throw(ArgumentError("When providing custom time-points for curriculum learning or \
            multiple shooting, each learning, each step/window i must include new
            measurements points. This does not hold for step $i. This means there are no \
            measurement points in the measurement table between \
            ($(tmax_prev), $(tmax_split)]."))
    end
    if maximum(splits[end]) < maximum(unique_t)
        throw(ArgumentError("When providing custom time-points for curriculum learning or \
            multiple shooting the last stage/window must cover all measurements in the
            measurement table. This does not hold as the end-point time for the last \
            interval is $(maximum(splits[end])) while the largest time-point in the \
            measurements table is $(maximum(unique_t))."))
    end
    if method == :multiple_shooting
        for i in 2:length(splits)
            splits[i][1] = maximum(splits[i - 1])
        end
        return splits
    end
    return _split_time_curriculum(splits, mdf, prob)
end

function _split_uniform_datapoints(prob::PEtabODEProblem, nsplits::Integer)
    mdf = prob.model_info.model.petab_tables[:measurements]
    mdf_sorted = mdf[sortperm(mdf.time), :]
    if nrow(mdf_sorted) < nsplits
        throw(ArgumentError("Number of data-points in the measurement table must be \
            greater than or equal to number of splits (number of curriculum steps). There \
            are $(nrow(mdf)) data-points and nsplits=$(nsplits)."))
    end

    out = Vector{Dict{Symbol, DataFrame}}(undef, nsplits)
    imaxs = _makechunks(collect(1:nrow(mdf)), nsplits) .|>
        maximum
    for i in 1:nsplits
        out[i] = copy(prob.model_info.model.petab_tables)
        out[i][:measurements] = mdf_sorted[1:imaxs[i], :]
    end
    return out
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
    splits = _makechunks(conditions_ids, nsplits)
    _split_condition_curriculum(splits, mdf, prob)
end

function _split_custom_conditions(prob::PEtabODEProblem,
        splits::Vector{<:Vector{T}}) where {T <: Union{String, Symbol}}
    splits_conditions = reduce(vcat, splits) .|> string
    conditions_ids = prob.model_info.simulation_info.conditionids[:simulation] .|> string
    for cid in conditions_ids
        cid in splits_conditions && continue
        throw(ArgumentError("The simulation condition $cid which has been defined for \
            the PEtab problem is not present in any of the custom condition curriculum \
            splits provided"))
    end
    mdf = prob.model_info.model.petab_tables[:measurements]
    _split_condition_curriculum(splits, mdf, prob)
end

function _get_unique_timepoints(mdf::DataFrame)::Vector{Float64}
    return mdf.time |>
           unique |>
           sort
end
