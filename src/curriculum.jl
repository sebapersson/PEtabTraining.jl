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
        model = __PEtabModel(original_problem.model_info.model.paths, petab_tables[i])
        petab_problems[i] = _PEtabODEProblem(model, original_problem)
    end
    return PEtabCurriculumProblem(split_algorithm, petab_problems, original_problem)
end
