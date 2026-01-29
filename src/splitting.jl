function _split_cl(split_alg::SplitTime, prob::PEtabODEProblem)
    _check_time_splits_cl(split_alg.spec, prob)

    unique_time_points = _get_unique_time_points(prob)
    if split_alg.spec isa Integer
        return last.(_makechunks(unique_time_points, split_alg.spec))
    end

    out = fill(0.0, length(split_alg.spec) + 1)
    for i in eachindex(split_alg.spec)
        idx_t = findall(t -> t < split_alg.spec[i], unique_time_points)
        out[i] = maximum(unique_time_points[idx_t])
    end
    out[end] = maximum(unique_time_points)
    return out
end
function _split_cl(split_alg::SplitData, prob::PEtabODEProblem)
    _check_data_splits_cl(split_alg.spec, prob)

    if split_alg.spec isa Integer
        measurements_df = prob.model_info.model.petab_tables[:measurements]
        return last.(_makechunks(collect(1:nrow(measurements_df)), split_alg.spec))
    end
    return split_alg.spec
end

function _split_uniform_time(prob::PEtabODEProblem, nsplits::Integer, method::Symbol)
    _check_n_splits(prob, nsplits, :time)

    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_t = _get_unique_time_points(prob.model_info.model.petab_tables[:measurements])
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
    unique_t = _get_unique_time_points(prob.model_info.model.petab_tables[:measurements])
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

function _split_uniform_datapoints(prob::PEtabODEProblem, nsplits::Integer)
    _check_n_splits(prob, nsplits, :datapoints)

    mdf_sorted = _get_sorted_measurements_df(prob)
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
    mdf_sorted = _get_sorted_measurements_df(prob)
    if splits[end] != nrow(mdf_sorted)
        throw(ArgumentError("When using custom splits based on data points, the number of \
            data points in the final interval must match the total number of data points \
            in the measurement table (expected: $(nrow(mdf_sorted)), provided \
            $(splits[end]))."))
    end

    return _split_curriculum(splits, mdf_sorted, prob, :datapoints)
end

function _check_n_splits(
        prob::PEtabODEProblem, nsplits::Integer, max_splits_criteria::Symbol)::Nothing
    if max_splits_criteria == :time
        mdf = prob.model_info.model.petab_tables[:measurements]
        max_splits = length(_get_unique_time_points(mdf))
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

function _check_time_splits_cl(n::Integer, prob::PEtabODEProblem)::Nothing
    unique_time_points = _get_unique_time_points(prob)
    if n > length(unique_time_points)
        throw(ArgumentError("Number of time splits $n exceeds number of unique time \
            points $(length(unique_time_points))."))
    end
    return nothing
end
function _check_time_splits_cl(split_points::Vector{<:Real}, prob::PEtabODEProblem)::Nothing
    if !issorted(split_points)
        throw(ArgumentError("Time split points must be sorted in ascending order."))
    end

    unique_time_points = _get_unique_time_points(prob)
    if maximum(unique_time_points) < split_points[end]
        throw(ArgumentError("The last time split point $(split_points[end]) is larger \
            than the maximum measurement time point $(maximum(unique_time_points))."))
    end

    # Check that each split will contain at least one time point
    unique_time_points = _get_unique_time_points(prob)
    t_start = unique_time_points[1]
    for (i, t_split) in pairs(split_points)
        if i != length(split_points)
            n_points = sum(
                (t_start .<= unique_time_points) .& (unique_time_points .< t_split)
            )
        else
            n_points = sum(
                (t_start .<= unique_time_points) .& (unique_time_points .< t_split)
            )
        end

        if n_points == 0
            interval = if i == length(split_points)
                "[($(t_start), $(t_split)]"
            else
                "[($(t_start), $(t_split))"
            end
            throw(ArgumentError("Time split point $t_split does not include any new \
                measurement time points, as there are no measurements between $(interval). \
                Each split must include at least one new measurement time point."))
        end
        t_start = t_split
    end
    return nothing
end

function _check_data_splits_cl(n::Integer, prob::PEtabODEProblem)::Nothing
    n_data_points = nrow(prob.model_info.model.petab_tables[:measurements])
    if n > n_data_points
        throw(ArgumentError("Number of data splits $n exceeds number of data points \
            $(n_data_points)."))
    end
    return nothing
end
function _check_data_splits_cl(
        split_points::Vector{<:Integer}, prob::PEtabODEProblem
    )::Nothing
    n_data_points = nrow(prob.model_info.model.petab_tables[:measurements])

    if !issorted(split_points)
        throw(ArgumentError("Data split points must be sorted in ascending order."))
    end


    if split_points[end] != n_data_points
        throw(ArgumentError("The last data split point $(split_points[end]) must be \
            equal to the number of data points $(n_data_points) in the measurement table."))
    end
    return nothing
end
