struct SplitUniform
    nsplits::Integer
    mode::Symbol
end
function SplitUniform(nsplits::Integer; mode::Symbol = :time)::SplitUniform
    @argcheck mode in [:time, :condition, :datapoints] "Mode must be either :time, \
         :datapoints or :condition"
    return SplitUniform(nsplits, mode)
end
struct SplitCustom{T}
    nsplits::Integer
    mode::Symbol
    splits::T
end
function SplitCustom(splits::Union{Vector{<:Real}, Vector{Vector{T}}};
        mode::Symbol = :time)::SplitCustom where {T <: Union{String, Symbol}}
    @argcheck mode in [:time, :condition, :datapoints] "Mode must be either :time, \
         :datapoints or :condition"

    # Internally, for a common interface, splits must be a Vector{Vector}.
    if mode == :condition &&
       !(splits isa Vector{Vector{Symbol}} || splits isa Vector{Vector{String}})
        throw(ArgumentError("If mode = :condition then provided splits must be a \
            Vector{Symbol} with condition to split over, e.g., [[:cond1, cond2], [:cond3]] \
            or a Vector{Vector{String}}, e.g, [[\"cond1\", \"cond2\"], [\"cond3\"]]))"))
    end
    if mode == :datapoints && !(splits isa Vector{<:Integer})
        throw(ArgumentError("If mode = :datapoints then provided splits must be a \
            Vector{Integer}, where entry i is number of data-points for curriculum step \
            i"))
    end
    if mode == :time && !(splits isa Vector{<:Real})
        throw(ArgumentError("If mode = :time then provided splits must be a \
            Vector{Real}, with the end-points for each curriculum step"))
    end

    for (i, split) in pairs(splits)
        !isempty(split) && continue
        throw(ArgumentError("Provided interval $i is empty. When providing custom \
            splits no interval must be empty"))
    end

    return SplitCustom(length(splits), mode, splits)
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
    mode = split_algorithm.mode
    if method == :multiple_shooting && (mode == :condition || mode == :datapoints)
        throw(ArgumentError("For multiple shotting, splitting windows over $mode \
            (mode = :$(mode)) is not allowed."))
    end

    if mode == :time
        return _split_custom_time(prob, split_algorithm.splits, method)
    elseif mode == :condition
        return _split_custom_conditions(prob, split_algorithm.splits)
    elseif mode == :datapoints
        return _split_custom_datapoints(prob, split_algorithm.splits)
    end
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

function _split_custom_time(prob::PEtabODEProblem, splits::Vector{<:Real}, method::Symbol)
    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_timepoints(mdf)
    n_time_points_included = 0
    for (i, tmax_split) in pairs(splits)
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
        splits = [[0.0, split] for split in splits]
        for i in 2:length(splits)
            splits[i][1] = maximum(splits[i - 1])
        end
        return splits
    end
    return _split_time_curriculum(splits, mdf, prob)
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
        splits::Vector{Vector{T}}) where {T <: Union{String, Symbol}}
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

# TODO: Common logic
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

function _split_custom_datapoints(prob::PEtabODEProblem, splits::Vector{<:Integer})
    nsplits = length(splits)
    for i in 1:(nsplits - 1)
        splits[i] < splits[i + 1] && continue
        throw(ArgumentError("When doing custom split based on data-points, each new \
            split interval must include more data-points than the previous. Interval \
            i = $(i+1) has less data-points than i = $i"))
    end

    mdf = prob.model_info.model.petab_tables[:measurements]
    mdf_sorted = mdf[sortperm(mdf.time), :]
    if splits[end] != nrow(mdf)
        throw(ArgumentError("When doing custom split based on data-points, the number \
            of data-points in the last interval must match the number of data-points in \
            the measurements table (in this case $(nrow(mdf)))"))
    end

    if nrow(mdf_sorted) < nsplits
        throw(ArgumentError("Number of data-points in the measurement table must be \
            greater than or equal to number of splits (number of curriculum steps). There \
            are $(nrow(mdf)) data-points and nsplits=$(nsplits)."))
    end

    out = Vector{Dict{Symbol, DataFrame}}(undef, nsplits)
    for i in 1:nsplits
        out[i] = copy(prob.model_info.model.petab_tables)
        out[i][:measurements] = mdf_sorted[1:splits[i], :]
    end
    return out
end

function _get_unique_timepoints(mdf::DataFrame)::Vector{Float64}
    return mdf.time |>
           unique |>
           sort
end
