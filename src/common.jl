
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

function _filter_condition_table!(petab_tables::Dict{Symbol, DataFrame})::Nothing
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
