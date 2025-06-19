using PEtab, PEtabTraining, Test

# Uniform splitting + curriculum
model_id = "Weber_BMC2015"
path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
@test_throws ArgumentError begin
    PEtabCurriculumProblem(petab_prob, SplitUniform(3; mode = :conditions))
end
@test_throws ArgumentError begin
    PEtabCurriculumProblem(petab_prob, SplitUniform(4; mode = :condition))
end
@test_throws ArgumentError begin
    PEtabCurriculumProblem(petab_prob, SplitUniform(14; mode = :time))
end
