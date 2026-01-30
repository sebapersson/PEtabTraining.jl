using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "helper.jl"))

function test_split_uniform_time(model_id, n_chunks)
    petab_prob = _get_petab_problem(model_id)
    cl_problem = PEtabClProblem(petab_prob, SplitTime(n_chunks))

    unique_t = PEtabTraining._get_unique_time_points(petab_prob)
    time_chunks = PEtabTraining._makechunks(unique_t, n_chunks)
    test_chunk_size(time_chunks)

    measurements_df = petab_prob.model_info.model.petab_tables[:measurements]
    for i in 1:n_chunks
        tmp_df = measurements_df[measurements_df[!, :time] .≤ maximum(time_chunks[i]), :]
        test_nllh(
            model_id, measurements_df, tmp_df, petab_prob, cl_problem.petab_problems[i]
        )
    end
    return nothing
end

function test_split_uniform_datapoints(model_id, n_chunks)
    petab_prob = _get_petab_problem(model_id)
    cl_problem = PEtabClProblem(petab_prob, SplitData(n_chunks))

    measurements_df = PEtabTraining._get_sorted_measurements_df(petab_prob)
    @test issorted(measurements_df.time)

    idx_chunk = last.(PEtabTraining._makechunks(collect(1:nrow(measurements_df)), n_chunks))
    for i in 1:n_chunks
        mdf_tmp = measurements_df[1:idx_chunk[i], :]
        test_nllh(
            model_id, measurements_df, mdf_tmp, petab_prob, cl_problem.petab_problems[i]
        )
    end
    return nothing
end

function test_split_custom_time(model_id, time_split)
    petab_prob = _get_petab_problem(model_id)
    cl_problem = PEtabClProblem(petab_prob, SplitTime(time_split))
    measurements_df = petab_prob.model_info.model.petab_tables[:measurements]

    @test length(cl_problem.petab_problems) == length(time_split) + 1

    for (i, max_val) in pairs(time_split)
        idx = findall(t -> t < max_val, measurements_df.time)
        mdf_tmp = measurements_df[idx, :]
        test_nllh(
            model_id, measurements_df, mdf_tmp, petab_prob, cl_problem.petab_problems[i]
        )
    end
    test_nllh(
        model_id, measurements_df, measurements_df, petab_prob, cl_problem.petab_problems[end]
    )
    return nothing
end

function test_split_custom_datapoints(model_id, chunk_sizes)
    petab_prob = _get_petab_problem(model_id)
    cl_problem = PEtabClProblem(petab_prob, SplitData(chunk_sizes))

    measurements_df = PEtabTraining._get_sorted_measurements_df(petab_prob)
    for (i, n_data_points) in pairs(chunk_sizes)
        mdf_tmp = measurements_df[1:n_data_points, :]
        test_nllh(
            model_id, measurements_df, mdf_tmp, petab_prob, cl_problem.petab_problems[i]
        )
    end
    return nothing
end

function test_output_regularization()
    petab_prob = _get_petab_problem("ude"; include_regularization = true)
    cl_prob = PEtabClProblem(petab_prob, SplitTime(4); regularization_obs = :reg_o)
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
        for n_chunks in [2, 3, 5]
            test_split_uniform_time("Boehm_JProteomeRes2014", n_chunks)
        end
        for n_chunks in [3, 4, 6]
            test_split_uniform_time("Weber_BMC2015", n_chunks)
        end
        for n_chunks in [3, 4]
            test_split_uniform_time("mm_julia", n_chunks)
        end
        for n_chunks in [2, 3]
            test_split_uniform_time("ude", n_chunks)
        end
        for n_chunks in [2, 7]
            test_split_uniform_time("Bachmann_MSB2011", n_chunks)
        end
    end

    @testset "Custom time" begin
        splits_test = [[15.0, 20.0, 100.0], [13.0, 105.0]]
        for time_split in splits_test
            test_split_custom_time("Boehm_JProteomeRes2014", time_split)
        end
        splits_test = [[1.5, 2.5], [1.3, 2.5, 4.0]]
        for time_split in splits_test
            test_split_custom_time("mm_julia", time_split)
        end
    end

    @testset "Uniform splitting data-points" begin
        for n_chunks in [2, 3, 5]
            test_split_uniform_datapoints("Boehm_JProteomeRes2014", n_chunks)
        end
        for n_chunks in [2, 3, 6]
            test_split_uniform_datapoints("mm_julia", n_chunks)
        end
        for n_chunks in [3, 4, 6]
            test_split_uniform_datapoints("Weber_BMC2015", n_chunks)
        end
        for n_chunks in [2, 7]
            test_split_uniform_datapoints("Bachmann_MSB2011", n_chunks)
        end
    end

    @testset "Custom data-points" begin
        splits_test = [[3, 6, 8, 48], [20, 21, 48]]
        for chunk_sizes in splits_test
            test_split_custom_datapoints("Boehm_JProteomeRes2014", chunk_sizes)
        end
        splits_test = [[3, 6, 22, 42], [20, 30, 42]]
        for chunk_sizes in splits_test
            test_split_custom_datapoints("mm_julia", chunk_sizes)
        end
    end

    @testset "Output regularization" begin
        test_output_regularization()
    end
end
