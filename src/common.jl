
function _makechunks(x::AbstractVector, n::Integer; overlap::Integer = 0)
    c = length(x) ÷ n
    chunks = [x[(1 + c * k):(k == n - 1 ? end : c * k + c)] for k in 0:(n - 1)]
    overlap == 0 && return chunks
    for i in 1:(length(chunks) - 1)
        chunks[i] = vcat(chunks[i], chunks[i + 1][1:overlap])
    end
    return chunks
end

function _get_specie_ids(prob::PEtabODEProblem)
    return _get_specie_ids(prob.model_info.model.speciemap)
end
function _get_specie_ids(speciemap)
    specie_ids = speciemap .|>
                 first .|>
                 string
    specie_ids = replace.(specie_ids, "(t)" => "")
    return specie_ids
end

function _filter_condition_table!(petab_tables::PEtab.PEtabTables)::Nothing
    conditions_df = petab_tables[:conditions]
    measurements_df = petab_tables[:measurements]
    cids_measurements_df = unique(measurements_df.simulationConditionId)
    cid_remove = String[]
    for cid in conditions_df.conditionId
        cid in cids_measurements_df && continue
        if "preequilibrationConditionId" in names(measurements_df)
            cond = cid in measurements_df.preequilibrationConditionId
            ismissing(cond) && continue
            cond == true && continue
        end
        push!(cid_remove, cid)
    end
    _conditions_df = filter(row -> !(row.conditionId in cid_remove), conditions_df)
    petab_tables[:conditions] = _conditions_df
    return nothing
end

function _get_unique_timepoints(mdf::DataFrame)::Vector{Float64}
    return mdf.time |>
           unique |>
           sort
end

function _get_measurements_df_sorted(prob::PEtabODEProblem)::DataFrame
    mdf = prob.model_info.model.petab_tables[:measurements]
    return mdf[sortperm(mdf.time), :]
end

function _splits_to_windows(splits::Vector{<:Real})
    splits = [[0.0, split] for split in splits]
    for i in 2:length(splits)
        splits[i][1] = maximum(splits[i - 1])
    end
    return splits
end

function _transform_x!(prob::PEtabODEProblem)::Nothing
    @unpack xnominal, xnominal_transformed, xnames, model_info = prob
    @views xnominal_transformed .= PEtab.transform_x(xnominal, xnames, model_info.xindices; to_xscale = true)
    return nothing
end

"""
_perm_from_labels(x::ComponentVector, y::ComponentVector)

Return a permutation ix such that getdata(y)[ix] == getdata(x).
"""
function _perm_from_labels(x::ComponentArrays.ComponentVector, y::ComponentArrays.ComponentVector)
    ix = fill(0, length(x))
    for (i, label) in pairs(ComponentArrays.labels(x))
        ix[i] = only(ComponentArrays.label2index(y, label))
    end
    return ix
end

function _check_regularization_obs(regularization_obs::Union{String, Symbol, Nothing},
        prob_original::PEtabODEProblem)::Nothing
    isnothing(regularization_obs) && return nothing
    regularization_obs = string(regularization_obs)
    measurements_original = prob_original.model_info.model.petab_tables[:measurements]
    observables_original = prob_original.model_info.model.petab_tables[:observables]
    @argcheck regularization_obs in observables_original.observableId "observableId \
        $(regularization_obs) must appear as an observable among the observables"
    @argcheck regularization_obs in measurements_original.observableId "observableId \
        $(regularization_obs) must appear for at least one measurement in the measurement \
        table"
    return nothing
end

function _check_regularization_specie(regularization_obs::Union{Nothing, String, Symbol},
        regularization_specie::Union{Nothing, String, Symbol})
    if !isnothing(regularization_specie) || !isnothing(regularization_obs)
        @argcheck !isnothing(regularization_obs) && !isnothing(regularization_specie) "If \
            regularization_obs is provided then regularization_specie must be provided"
        regularization_obs = string(regularization_obs)
        regularization_specie = string(regularization_specie)
    end
    return nothing
end

_string(x::String) = x
_string(::Nothing) = nothing
_string(x::Symbol) = string(x)
