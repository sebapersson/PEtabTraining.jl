struct PEtabCurriculumProblem
    split_algorithm::Any
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
end
function PEtabCurriculumProblem(
        prob_original::PEtabODEProblem, split_algorithm)::PEtabCurriculumProblem
    petab_tables = _split(split_algorithm, prob_original, :curriculum)
    petab_problems = Vector{PEtabODEProblem}(undef, split_algorithm.nsplits)
    for i in eachindex(petab_tables)
        _filter_condition_table!(petab_tables[i])
        model = PEtab._PEtabModel(prob_original.model_info.model.paths,
            petab_tables[i], false, false, true, false)
        petab_problems[i] = _PEtabODEProblem(model, prob_original)
    end
    return PEtabCurriculumProblem(split_algorithm, petab_problems, prob_original)
end
function _split_curriculum(splits, mdf::DataFrame, prob::PEtabODEProblem,
        mode::Symbol)::Vector{Dict{Symbol, DataFrame}}
    @assert mode in [:datapoints, :time, :conditions]
    out = Vector{Dict{Symbol, DataFrame}}(undef, length(splits))
    if mode in [:datapoints, :time]
        for (i, split) in pairs(splits)
            out[i] = copy(prob.model_info.model.petab_tables)
            if mode == :time
                tmax = maximum(split)
                out[i][:measurements] = mdf[mdf[!, :time] .≤ tmax, :]
            else
                out[i][:measurements] = mdf[1:maximum(splits[i]), :]
            end
        end
    elseif mode == :conditions
        for (i, split) in pairs(splits)
            out[i] = copy(prob.model_info.model.petab_tables)
            irow = findall(x -> x in split, mdf.simulationConditionId)
            out[i][:measurements] = mdf[irow, :]
        end
    end
    return out
end
