using CSV, DataFrames, PEtab, PEtabTraining, Test

function test_split_uniform_time(model_id, n_stages)
    path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
    petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
    stage_problems = PEtabCurriculumProblem(petab_prob, SplitUniform(n_stages))

    # Check that the chunking produces expected chunk-sizes
    mdf = petab_prob.model_info.model.petab_tables[:measurements]
    unique_t = mdf.time |> sort |> unique
    time_chunks = PEtabTraining._makechunks(unique_t, n_stages)
    for i in 1:(length(time_chunks) - 1)
        # Sometimes chunks cannot be divided uniformly, then the last is bigger
        if i == length(time_chunks) - 1
            @test abs(length(time_chunks[i]) - length(time_chunks[i + 1])) < n_stages
        else
            @test length(time_chunks[i]) == length(time_chunks[i + 1])
        end
    end

    # Check that the likelihood is correct
    for i in 1:n_stages
        mdf_tmp = mdf[mdf[!, :time] .≤ maximum(time_chunks[i]), :]
        CSV.write(petab_prob.model_info.model.paths[:measurements], mdf_tmp, delim = '\t')
        petab_prob_ref = PEtabModel(path_yaml) |> PEtabODEProblem
        CSV.write(petab_prob.model_info.model.paths[:measurements], mdf, delim = '\t')
        nllh_ref = petab_prob_ref.nllh(get_x(petab_prob_ref))
        nllh_test = stage_problems.petab_problems[i].nllh(get_x(petab_prob))
        @test nllh_ref == nllh_test
    end
    return nothing
end

function test_split_uniform_conditions(model_id, n_stages)
    path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
    petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
    stage_problems = PEtabCurriculumProblem(
        petab_prob, SplitUniform(n_stages; mode = :condition))

    # Check that the chunking produces expected chunk-sizes
    mdf = petab_prob.model_info.model.petab_tables[:measurements]
    unique_conditions = mdf.simulationConditionId |> unique
    condition_chunks = PEtabTraining._makechunks(unique_conditions, n_stages)
    for i in 1:(length(condition_chunks) - 1)
        # Sometimes chunks cannot be divided uniformly, then the last is bigger
        if i == length(condition_chunks) - 1
            @test abs(length(condition_chunks[i]) - length(condition_chunks[i + 1])) <
                  n_stages
        else
            @test length(condition_chunks[i]) == length(condition_chunks[i + 1])
        end
    end

    # Check that the likelihood is correct
    for i in 1:n_stages
        irow = findall(x -> x in condition_chunks[i], mdf.simulationConditionId)
        mdf_tmp = mdf[irow, :]
        CSV.write(petab_prob.model_info.model.paths[:measurements], mdf_tmp, delim = '\t')
        petab_prob_ref = PEtabModel(path_yaml) |> PEtabODEProblem
        CSV.write(petab_prob.model_info.model.paths[:measurements], mdf, delim = '\t')
        nllh_ref = petab_prob_ref.nllh(get_x(petab_prob_ref))
        nllh_test = stage_problems.petab_problems[i].nllh(get_x(petab_prob))
        @test nllh_ref == nllh_test
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

@testset "Uniform splitting conditions" begin
    for n_stages in [2]
        test_split_uniform_conditions("Weber_BMC2015", n_stages)
    end
    for n_stages in [2, 5, 7]
        test_split_uniform_conditions("Bachmann_MSB2011", n_stages)
    end
end
