struct PEtabCurriculumProblem
    split_algorithm::Any
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
end
function PEtabCurriculumProblem(
        original_problem::PEtabODEProblem, split_algorithm)::PEtabCurriculumProblem
    petab_tables = _split(split_algorithm, original_problem)
    petab_problems = Vector{PEtabODEProblem}(undef, split_algorithm.nsplits)
    for i in eachindex(petab_tables)
        model = PEtab._PEtabModel(original_problem.model_info.model.paths,
            petab_tables[i], false, false, true, false)
        petab_problems[i] = _PEtabODEProblem(model, original_problem)
    end
    return PEtabCurriculumProblem(split_algorithm, petab_problems, original_problem)
end
function _split_time_curriculum(
        splits, mdf::DataFrame, prob::PEtabODEProblem)::Vector{Dict{Symbol, DataFrame}}
    out = Vector{Dict{Symbol, DataFrame}}(undef, length(splits))
    for (i, split) in pairs(splits)
        tmax = maximum(split)
        out[i] = copy(prob.model_info.model.petab_tables)
        out[i][:measurements] = mdf[mdf[!, :time] .≤ tmax, :]
    end
    return out
end

function _split_condition_curriculum(
        splits, mdf::DataFrame, prob::PEtabODEProblem)::Vector{Dict{Symbol, DataFrame}}
    out = Vector{Dict{Symbol, DataFrame}}(undef, length(splits))
    for (i, split) in pairs(splits)
        out[i] = copy(prob.model_info.model.petab_tables)
        irow = findall(x -> x in split, mdf.simulationConditionId)
        out[i][:measurements] = mdf[irow, :]
    end
    return out
end
