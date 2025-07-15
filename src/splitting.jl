struct SplitUniform
    nsplits::Integer
    mode::Symbol
end
function SplitUniform(nsplits::Integer; mode::Symbol = :time)::SplitUniform
    _check_mode(mode)
    return SplitUniform(nsplits, mode)
end
struct SplitCustom{T}
    nsplits::Integer
    mode::Symbol
    splits::T
end
function SplitCustom(splits::Union{Vector{<:Real}, Vector{Vector{T}}};
        mode::Symbol = :time)::SplitCustom where {T <: Union{String, Symbol}}
    _check_mode(mode)

    # Sanity check correct user input
    s = splits
    if mode == :time && !(s isa Vector{<:Real})
        throw(ArgumentError("If `mode = :time`, then `splits` must be a `Vector{<:Real}`, \
            where each entry represents the end-point of a curriculum step or multiple \
            shooting window, e.g., [4.0, 20.0, 30.0]."))
    end
    if mode == :condition && !(s isa Vector{Vector{Symbol}} || s isa Vector{Vector{String}})
        throw(ArgumentError("If `mode = :condition`, then `splits` must be a \
             `Vector{Vector}` of simulation condition groups to split over, \
             e.g., [[:cond1, :cond2], [:cond3]]."))
    end
    if mode == :datapoints && !(s isa Vector{<:Integer})
        throw(ArgumentError("If `mode = :datapoints`, then splits must be a \
            `Vector{Integer}`, where each entry specifies the number of data points from
            the measurement table to include in each curriculum step, e.g., [5, 10, 20]"))
    end
    for (i, split) in pairs(splits)
        !isempty(split) && continue
        throw(ArgumentError("Interval $i is empty in the provided splits. When \
            providing custom splits, no interval may be empty."))
    end

    return SplitCustom(length(splits), mode, splits)
end

function _split(split_algorithm::SplitUniform, prob::PEtabODEProblem, method::Symbol)
    _check_mode(split_algorithm.mode, method)

    if split_algorithm.mode == :time
        return _split_uniform_time(prob, split_algorithm.nsplits, method)
    elseif split_algorithm.mode == :condition
        return _split_uniform_conditions(prob, split_algorithm.nsplits)
    elseif split_algorithm.mode == :datapoints
        return _split_uniform_datapoints(prob, split_algorithm.nsplits)
    end
end
function _split(split_algorithm::SplitCustom, prob::PEtabODEProblem, method::Symbol)
    _check_mode(split_algorithm.mode, method)

    if split_algorithm.mode == :time
        return _split_custom_time(prob, split_algorithm.splits, method)
    elseif split_algorithm.mode == :condition
        return _split_custom_conditions(prob, split_algorithm.splits)
    elseif split_algorithm.mode == :datapoints
        return _split_custom_datapoints(prob, split_algorithm.splits)
    end
end

function _split_uniform_time(prob::PEtabODEProblem, nsplits::Integer, method::Symbol)
    _check_n_splits(prob, nsplits, :time)

    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_timepoints(prob.model_info.model.petab_tables[:measurements])
    if method == :multiple_shooting
        return _makechunks(unique_t, nsplits; overlap = 1)
    end
    splits = _makechunks(unique_t, nsplits; overlap = 0)
    return _split_curriculum(splits, mdf, prob, :time)
end

function _split_custom_time(prob::PEtabODEProblem, splits::Vector{<:Real}, method::Symbol)
    _check_n_splits(prob, length(splits), :time)

    # Sanity check user provided reasonable splitting intervals
    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_timepoints(prob.model_info.model.petab_tables[:measurements])
    n_time_points_included = 0
    for (i, tmax_split) in pairs(splits)
        n_time_points_split = sum(unique_t .≤ tmax_split)
        if n_time_points_included < n_time_points_split
            n_time_points_included = n_time_points_split
            continue
        end

        tmax_prev = i == 1 ? 0.0 : maximum(splits[i - 1])
        throw(ArgumentError("Invalid custom time split: split $i does not include any \
            new measurement points. No measurement times are found in the interval \
            ($(tmax_prev), $(tmax_split)]. Each curriculum step or multiple shooting \
             window must include at least one new time point."))
    end
    if maximum(splits[end]) < maximum(unique_t)
        throw(ArgumentError("Invalid custom time split: the final curriculum step or \
            window must cover all measurement time points. The last interval ends at \
            $(maximum(splits[end])), but the latest time point in the measurement table \
            is $(maximum(unique_t)). Provide splits such that the last time points is \
            ≥ $(maximum(unique_t))"))
    end

    if method == :multiple_shooting
        return _splits_to_windows(splits)
    end
    return _split_curriculum(splits, mdf, prob, :time)
end

function _split_uniform_conditions(prob::PEtabODEProblem, nsplits::Integer)
    _check_n_splits(prob, nsplits, :conditions)

    conditions_ids = string.(prob.model_info.simulation_info.conditionids[:simulation])
    mdf = prob.model_info.model.petab_tables[:measurements]
    splits = _makechunks(conditions_ids, nsplits)
    _split_curriculum(splits, mdf, prob, :conditions)
end

function _split_custom_conditions(prob::PEtabODEProblem,
        splits::Vector{Vector{T}}) where {T <: Union{String, Symbol}}
    _check_n_splits(prob, length(splits), :conditions)
    splits_conditions = string.(reduce(vcat, splits))
    conditions_ids = string.(prob.model_info.simulation_info.conditionids[:simulation])
    for cid in conditions_ids
        cid in splits_conditions && continue
        throw(ArgumentError("The simulation condition $cid, which is defined in the \
            PEtab problem, is not present in any of the custom condition curriculum \
            splits provided."))
    end

    mdf = prob.model_info.model.petab_tables[:measurements]
    _split_curriculum(splits, mdf, prob, :conditions)
end

function _split_uniform_datapoints(prob::PEtabODEProblem, nsplits::Integer)
    _check_n_splits(prob, nsplits, :datapoints)

    mdf_sorted = _get_measurements_df_sorted(prob)
    splits = _makechunks(collect(1:nrow(mdf_sorted)), nsplits)
    return _split_curriculum(splits, mdf_sorted, prob, :datapoints)
end

function _split_custom_datapoints(prob::PEtabODEProblem, splits::Vector{<:Integer})
    nsplits = length(splits)
    _check_n_splits(prob, nsplits, :datapoints)
    for i in 1:(nsplits - 1)
        splits[i] < splits[i + 1] && continue
        throw(ArgumentError("When using custom splits based on data points, each new \
            interval must include more data points than the previous one. Interval \
            $(i + 1) has fewer ($(splits[i+1]) data points) data points than \
             $i ($(splits[i]) data points)."))
    end
    mdf_sorted = _get_measurements_df_sorted(prob)
    if splits[end] != nrow(mdf_sorted)
        throw(ArgumentError("When using custom splits based on data points, the number of \
            data points in the final interval must match the total number of data points \
            in the measurement table (expected: $(nrow(mdf_sorted)), provided \
            $(splits[end]))."))
    end

    return _split_curriculum(splits, mdf_sorted, prob, :datapoints)
end

function _check_mode(mode::Symbol)::Nothing
    @argcheck mode in [:time, :datapoints, :condition] "Splitting mode must be one of \
        :time, :datapoints, or :condition."
    return nothing
end
function _check_mode(mode::Symbol, method::Symbol)::Nothing
    if method == :multiple_shooting && (mode == :condition || mode == :datapoints)
        throw(ArgumentError("For multiple shooting, splitting over `$mode` \
            (mode = :$(mode)) is not allowed."))
    end
    return nothing
end

function _check_n_splits(
        prob::PEtabODEProblem, nsplits::Integer, max_splits_criteria::Symbol)::Nothing
    if max_splits_criteria == :time
        mdf = prob.model_info.model.petab_tables[:measurements]
        max_splits = length(_get_unique_timepoints(mdf))
        str1 = "unique time-points"
    elseif max_splits_criteria == :datapoints
        mdf = prob.model_info.model.petab_tables[:measurements]
        max_splits = nrow(mdf)
        str1 = "data-points in the measurement table"
    elseif max_splits_criteria == :conditions
        simulation_ids = prob.model_info.simulation_info.conditionids[:simulation]
        max_splits = length(simulation_ids)
        str1 = "simulation conditions"
    end
    if nsplits ≤ max_splits
        return nothing
    end
    throw(ArgumentError("The number of $str1 must be smaller than or equal to the number \
        of splits (curriculum stages or multiple shooting windows). There are \
        $(max_splits) $str1, but nsplits = $(nsplits)."))
end
