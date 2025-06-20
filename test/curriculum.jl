using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "helper.jl"))

function test_split_uniform_time(model_id, n_stages)
    path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
    petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
    stage_problems = PEtabCurriculumProblem(petab_prob, SplitUniform(n_stages))

    mdf = petab_prob.model_info.model.petab_tables[:measurements]
    unique_t = PEtabTraining._get_unique_timepoints(mdf)
    time_chunks = PEtabTraining._makechunks(unique_t, n_stages)
    test_chunk_size(time_chunks)

    for i in 1:n_stages
        mdf_tmp = mdf[mdf[!, :time] .≤ maximum(time_chunks[i]), :]
        test_nllh(path_yaml, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_split_custom_time(model_id, splits)
    path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
    petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
    stage_problems = PEtabCurriculumProblem(petab_prob, SplitCustom(splits))
    mdf = petab_prob.model_info.model.petab_tables[:measurements]

    for (i, max_val) in pairs(splits)
        mdf_tmp = mdf[mdf[!, :time] .≤ max_val, :]
        test_nllh(path_yaml, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_split_uniform_conditions(model_id, n_stages)
    path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
    petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
    stage_problems = PEtabCurriculumProblem(petab_prob, SplitUniform(n_stages; mode = :condition))

    mdf = petab_prob.model_info.model.petab_tables[:measurements]
    unique_conditions = mdf.simulationConditionId |> unique
    condition_chunks = PEtabTraining._makechunks(unique_conditions, n_stages)
    test_chunk_size(condition_chunks)

    for i in 1:n_stages
        irow = findall(x -> x in condition_chunks[i], mdf.simulationConditionId)
        mdf_tmp = mdf[irow, :]
        test_nllh(path_yaml, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

@testset "Uniform splitting time" begin
    for n_stages in [2, 3, 5]
        test_split_uniform_time("Boehm_JProteomeRes2014", n_stages)
    end
    for n_stages in [3, 4, 6]
        test_split_uniform_time("Weber_BMC2015", n_stages)
    end
    for n_stages in [2, 7]
        test_split_uniform_time("Bachmann_MSB2011", n_stages)
    end
end

@testset "Uniform custom time" begin
    splits_test = [[15.0, 20.0, 100.0, 240.0], [13.0, 25.0, 105.0, 250.0]]
    for splits in splits_test
        test_split_custom_time("Boehm_JProteomeRes2014", splits)
    end
end

@testset "Uniform splitting conditions" begin
    for n_stages in [2]
        test_split_uniform_conditions("Weber_BMC2015", n_stages)
    end
    for n_stages in [2, 5, 7]
        test_split_uniform_conditions("Bachmann_MSB2011", n_stages)
    end
end

model_id = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem

splits =
stage_problems = PEtabCurriculumProblem(petab_prob, SplitCustom(splits))
