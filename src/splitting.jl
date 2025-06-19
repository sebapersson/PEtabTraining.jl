struct SplitUniform
    nsplits::Integer
    mode::Symbol
end
function SplitUniform(nsplits::Integer; mode::Symbol = :time)::SplitUniform
    @argcheck mode in [:time, :condition]
    return SplitUniform(nsplits, mode)
end

function _split(split_algorithm::SplitUniform,
        prob::PEtabODEProblem)::Vector{Dict{Symbol, DataFrame}}
    if split_algorithm.mode == :time
        return _split_uniform_time(prob, split_algorithm.nsplits)
    end
    return _split_uniform_conditions(prob, split_algorithm.nsplits)
end

function _split_uniform_time(
        prob::PEtabODEProblem, nsplits::Integer)::Vector{Dict{Symbol, DataFrame}}
    mdf = prob.model_info.model.petab_tables[:measurements]
    unique_timepoints = mdf.time |> unique |> sort
    if length(unique_timepoints) < nsplits
        throw(ArgumentError("Number of unique time-points in the measurement table must be \
          greater than or equal to number of splits (number of curriculum stages or \
          multiple shooting windows). There are $(length(unique_timepoints)) unique \
          time-points and nsplits=$(nsplits)."),
        )
    end

    tmaxs = _makechunks(unique_timepoints, nsplits) .|> maximum
    out = Vector{Dict{Symbol, DataFrame}}(undef, nsplits)
    for i in 1:nsplits
        out[i] = copy(prob.model_info.model.petab_tables)
        out[i][:measurements] = mdf[mdf[!, :time] .≤ tmaxs[i], :]
    end
    return out
end

function _split_uniform_conditions(
        prob::PEtabODEProblem, nsplits::Integer)::Vector{Dict{Symbol, DataFrame}}
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

struct SplitCustom{T <: Union{AbstractFloat, String, Symbol}}
    nsplits::Integer
    mode::Symbol
    splits::Vector{Vector{T}}
end
