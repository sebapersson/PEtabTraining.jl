struct SplitUniform
    nsplits::Integer
    mode::Symbol
end
function SplitUniform(nsplits::Integer; mode::Symbol = :time)::SplitUniform
    @argcheck mode in [:time, :condition]
    return SplitUniform(nsplits, mode)
end

function _split(split_algorithm::SplitUniform, prob::PEtabODEProblem)::Vector{Dict{Symbol, DataFrame}}
    if split_algorithm.mode == :time
        return _split_uniform_time(prob, split_algorithm.nsplits)
    end
end

function _split_uniform_time(prob::PEtabODEProblem, nsplits::Integer)::Vector{Dict{Symbol, DataFrame}}
    measurements_df = prob.model_info.model.petab_tables[:measurements]
    unique_timepoints = measurements_df.time |>
           unique |>
           sort
    if length(unique_timepoints) < nsplits
        ArgumentError("Number of unique time-points in the measurement table must be \
            greater than or equal to number of splits (number of curriculum stages or \
            multiple shooting windows). There are $(length(unique_timepoints)) unique \
            time-points and nsplits=$(nsplits).")
    end

    tmaxs = _makechunks(unique_timepoints, nsplits) .|>
            maximum
    out = Vector{Dict{Symbol, DataFrame}}(undef, nsplits)
    for i in 1:nsplits
        out[i] = copy(prob.model_info.model.petab_tables)
        out[i][:measurements] = measurements_df[measurements_df[!, :time] .≤ tmaxs[i], :]
    end
    return out
end

struct SplitCustom{T <: Union{AbstractFloat, String, Symbol}}
    nsplits::Integer
    mode::Symbol
    splits::Vector{Vector{T}}
end
function SplitCustom(splits::Vector{Vector{T}}; mode::Symbol = :time)::SplitCustom where T <: Union{Real, String, Symbol}
    @argcheck mode in [:time, :condition]
    if mode == :time
        @argcheck eltype(splits) isa Vector{<:Real}
    else
        @argcheck all(eltype.(splits) .== Symbol) || all(eltype.(splits) .== String)
    end
    return SplitUniform(length(splits), mode, splits)
end
