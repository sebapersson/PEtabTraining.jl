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

function _split_ms(split_alg::SplitTime, prob::PEtabODEProblem)
    _check_time_splits_ms(split_alg.spec, prob)

    unique_time_points = _get_unique_time_points(prob)
    if split_alg.spec isa Integer
        windows = _makechunks(unique_time_points, split_alg.spec; overlap = 1)
        return [[first(window), last(window)] for window in windows]
    end

    window_splits = [unique_time_points[1], split_alg.spec..., unique_time_points[end]]
    return [[window_splits[i], window_splits[i + 1]] for i in 1:(length(window_splits) - 1)]
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
        throw(ArgumentError("Time split points must be sorted in ascending order, does \
            not hold for $(split_points)"))
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
            interval = "[$(t_start), $(t_split))"
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
        throw(ArgumentError("Data split points must be sorted in ascending order, does \
            not hold for $(split_points)."))
    end


    if split_points[end] != n_data_points
        throw(ArgumentError("The last data split point $(split_points[end]) must be \
            equal to the number of data points $(n_data_points) in the measurement table."))
    end
    return nothing
end

function _check_time_splits_ms(n::Integer, prob::PEtabODEProblem)::Nothing
    return _check_data_splits_cl(n, prob)
end
function _check_time_splits_ms(split_points::Vector{<:Real}, prob::PEtabODEProblem)::Nothing
    if !issorted(split_points)
        throw(ArgumentError("Time split points must be sorted in ascending order, does \
            not hold for $(split_points)"))
    end

    unique_time_points = _get_unique_time_points(prob)
    if maximum(unique_time_points) ≤ split_points[end]
        throw(ArgumentError("The last time split point $(split_points[end]) is larger or \
            equal to the maximum measurement time point $(maximum(unique_time_points))."))
    end

    if minimum(unique_time_points) ≥ split_points[1]
        throw(ArgumentError("The first time split point $(split_points[1]) is smaller or \
            equal to the minimum measurement time point $(minimum(unique_time_points))."))
    end
    return nothing
end
