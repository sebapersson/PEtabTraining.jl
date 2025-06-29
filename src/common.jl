
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
