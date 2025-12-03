using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "helper.jl"))

function test_split_uniform_time(model_id, n_stages)
    petab_prob = _get_petab_problem(model_id)
    stage_problems = PEtabCLProblem(petab_prob, SplitUniform(n_stages))

    mdf = petab_prob.model_info.model.petab_tables[:measurements]
    unique_t = PEtabTraining._get_unique_timepoints(mdf)
    time_chunks = PEtabTraining._makechunks(unique_t, n_stages)
    test_chunk_size(time_chunks)

    for i in 1:n_stages
        mdf_tmp = mdf[mdf[!, :time] .≤ maximum(time_chunks[i]), :]
        test_nllh(model_id, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_split_uniform_datapoints(model_id, n_stages)
    petab_prob = _get_petab_problem(model_id)
    stage_problems = PEtabCLProblem(
        petab_prob, SplitUniform(n_stages; mode = :datapoints))

    mdf_sorted = PEtabTraining._get_measurements_df_sorted(petab_prob)
    @test issorted(mdf_sorted.time)
    imaxs = PEtabTraining._makechunks(collect(1:nrow(mdf_sorted)), n_stages) .|>
            maximum
    for i in 1:n_stages
        mdf_tmp = mdf_sorted[1:imaxs[i], :]
        test_nllh(
            model_id, mdf_sorted, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_split_custom_time(model_id, splits)
    petab_prob = _get_petab_problem(model_id)
    stage_problems = PEtabCLProblem(petab_prob, SplitCustom(splits; mode = :time))
    mdf = petab_prob.model_info.model.petab_tables[:measurements]

    for (i, max_val) in pairs(splits)
        mdf_tmp = mdf[mdf[!, :time] .≤ max_val, :]
        test_nllh(model_id, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_split_uniform_conditions(model_id, n_stages)
    petab_prob = _get_petab_problem(model_id)
    stage_problems = PEtabCLProblem(
        petab_prob, SplitUniform(n_stages; mode = :condition))

    mdf = petab_prob.model_info.model.petab_tables[:measurements]
    unique_conditions = mdf.simulationConditionId |> unique
    condition_chunks = PEtabTraining._makechunks(unique_conditions, n_stages)
    test_chunk_size(condition_chunks)

    for i in 1:n_stages
        irow = findall(x -> x in condition_chunks[i], mdf.simulationConditionId)
        mdf_tmp = mdf[irow, :]
        test_nllh(model_id, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_split_conditions_custom(model_id, splits)
    petab_prob = _get_petab_problem(model_id)
    stage_problems = PEtabCLProblem(
        petab_prob, SplitCustom(splits; mode = :condition))
    mdf = petab_prob.model_info.model.petab_tables[:measurements]

    for (i, split) in pairs(splits)
        irow = findall(x -> x in split, mdf.simulationConditionId)
        mdf_tmp = mdf[irow, :]
        test_nllh(model_id, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_split_custom_datapoints(model_id, splits)
    petab_prob = _get_petab_problem(model_id)
    stage_problems = PEtabCLProblem(
        petab_prob, SplitCustom(splits; mode = :datapoints))
    mdf = petab_prob.model_info.model.petab_tables[:measurements]
    mdf_sorted = mdf[sortperm(mdf.time), :]
    for (i, imax) in pairs(splits)
        mdf_tmp = mdf_sorted[1:imax, :]
        test_nllh(model_id, mdf, mdf_tmp, petab_prob, stage_problems.petab_problems[i])
    end
    return nothing
end

function test_output_regularization()
    petab_prob = _get_petab_problem("ude")
    cl_prob = PEtabCLProblem(petab_prob, SplitUniform(4); regularization_obs = :reg_o)
    for prob in cl_prob.petab_problems
        measurements_stage = prob.model_info.model.petab_tables[:measurements]
        @test "reg_o" in measurements_stage.observableId
        i_rows = findall(x -> x == "reg_o", measurements_stage.observableId)
        @test length(i_rows) == 1
        @test measurements_stage.time[i_rows[1]] == maximum(measurements_stage.time)
    end
    return nothing
end

@testset "Curriculum learning" begin
    @testset "Uniform splitting time" begin
        for n_stages in [2, 3, 5]
            test_split_uniform_time("Boehm_JProteomeRes2014", n_stages)
        end
        for n_stages in [3, 4, 6]
            test_split_uniform_time("Weber_BMC2015", n_stages)
        end
        for n_stages in [3, 4]
            test_split_uniform_time("mm_julia", n_stages)
        end
        for n_stages in [2, 3]
            test_split_uniform_time("ude", n_stages)
        end
        for n_stages in [2, 7]
            test_split_uniform_time("Bachmann_MSB2011", n_stages)
        end
    end

    @testset "Custom time" begin
        splits_test = [[15.0, 20.0, 100.0, 240.0], [13.0, 25.0, 105.0, 250.0]]
        for splits in splits_test
            test_split_custom_time("Boehm_JProteomeRes2014", splits)
        end
        splits_test = [[1.5, 2.5, 10.0], [1.3, 2.5, 4.0, 10.0]]
        for splits in splits_test
            test_split_custom_time("mm_julia", splits)
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

    @testset "Custom splitting conditions" begin
        model_id = "Fujita_SciSignal2010"
        path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
        petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
        splits = [[:condition_step_00_1, :condition_step_00_3],
            [:condition_step_01_0, :condition_step_03_0],
            [:condition_step_10_0, :condition_step_30_0]]
        splits_str = [string.(split) for split in splits]
        test_split_conditions_custom(model_id, splits)
        test_split_conditions_custom(model_id, splits_str)
    end

    @testset "Uniform splitting data-points" begin
        for n_stages in [2, 3, 5]
            test_split_uniform_datapoints("Boehm_JProteomeRes2014", n_stages)
        end
        for n_stages in [2, 3, 6]
            test_split_uniform_datapoints("Boehm_JProteomeRes2014", n_stages)
        end
        for n_stages in [3, 4, 6]
            test_split_uniform_datapoints("Weber_BMC2015", n_stages)
        end
        for n_stages in [2, 7]
            test_split_uniform_datapoints("Bachmann_MSB2011", n_stages)
        end
    end

    @testset "Custom data-points" begin
        splits_test = [[3, 6, 8, 48], [20, 21, 30, 48]]
        for splits in splits_test
            test_split_custom_datapoints("Boehm_JProteomeRes2014", splits)
        end
        splits_test = [[3, 6, 22, 42], [20, 30, 42]]
        for splits in splits_test
            test_split_custom_datapoints("mm_julia", splits)
        end
    end

    @testset "Output regularization" begin
        test_output_regularization()
    end
end
