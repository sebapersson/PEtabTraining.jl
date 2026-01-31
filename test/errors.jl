using PEtab, PEtabTraining, Test

# Uniform splitting + curriculum
model_id = "Weber_BMC2015"
path_yaml = joinpath(@__DIR__, "models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
@test_throws ArgumentError begin
    PEtabClProblem(petab_prob, SplitTime(14))
end
@test_throws ArgumentError begin
    PEtabClProblem(petab_prob, SplitData(137))
end

# Custom splitting + curriculum
model_id = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem

@test_throws ArgumentError begin
    PEtabClProblem(petab_prob, SplitTime([3.0, 2.0, 4.0]))
end
@test_throws ArgumentError begin
    PEtabClProblem(petab_prob, SplitTime([3.0, 4.0, 250]))
end
@test_throws ArgumentError begin
    PEtabClProblem(petab_prob, SplitTime([3.0, 12.0, 13.0, 230]))
end
@test_throws ArgumentError begin
    PEtabClProblem(petab_prob, SplitData([3, 4, 49]))
end
@test_throws ArgumentError begin
    PEtabClProblem(petab_prob, SplitData([3, 2, 48]))
end

# Multiple shooting unique errors
model_id = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
@test_throws ArgumentError begin
    PEtabMsProblem(petab_prob, SplitTime(49))
end
@test_throws ArgumentError begin
    PEtabMsProblem(petab_prob, SplitTime([3.0, 240.0]))
end
@test_throws ArgumentError begin
    PEtabMsProblem(petab_prob, SplitTime([-3.0, 230.0]))
end

model_id = "Weber_BMC2015"
path_yaml = joinpath(@__DIR__, "models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
@test_throws ArgumentError begin
    PEtabMsProblem(petab_prob, SplitTime(4))
end

# Combined approach
model_id = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
prob_cl_ms = PEtabClMsProblem(petab_prob, SplitTime(3))
x_test = get_x(prob_cl_ms.petab_problems[1])
@test_throws ArgumentError begin
    set_u0_ms_windows!(x_test, prob_cl_ms, 4; init = MsInitConstant(5.0))
end
