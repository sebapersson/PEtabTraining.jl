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

# Custom splitting + curriculum
model_id = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
@test_throws ArgumentError begin
    PEtabCurriculumProblem(petab_prob, SplitCustom([3.0, 2.0, 4.0]))
end
@test_throws ArgumentError begin
    PEtabCurriculumProblem(petab_prob, SplitCustom([3.0, 4.0, 250]))
end
@test_throws ArgumentError begin
    PEtabCurriculumProblem(petab_prob, SplitCustom([3.0, 10.0, 230]))
end
model_id = "Fujita_SciSignal2010"
path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
@test_throws ArgumentError begin
    splits = [[:condition_step_00_1, :condition_step_00_3], [:condition_step_01_0, :condition_step_03_0], [:condition_step_30_0]]
    PEtabCurriculumProblem(petab_prob, SplitCustom(splits))
end
@test_throws ArgumentError begin
    splits = [[:condition_step_00_1, :condition_step_00_3], [:condition_step_01_0, :condition_step_03_0], Symbol[]]
    SplitCustom(splits)
end
