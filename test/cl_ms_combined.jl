using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "common.jl"))

function test_cl_ms(model_id, split_alg::SplitTime; window_u0_scale = :lin)
    @assert window_u0_scale in [:lin, :log10]

    prob_original = _get_petab_problem(model_id)
    prob_cl_ms = PEtabClMsProblem(
        prob_original, split_alg; window_u0_scale = window_u0_scale
    )

    # Test that windows are created correctly
    for i in length(prob_cl_ms.ms_windows):-1:2
        windows_stage = prob_cl_ms.ms_windows[Symbol("stage$i")]
        windows_stage_prev = prob_cl_ms.ms_windows[Symbol("stage$(i - 1)")]
        for j in eachindex(windows_stage)
            test_value = unique(reduce(vcat, windows_stage_prev[j:(j + 1)]))
            test_value = [first(test_value), last(test_value)]
            @test windows_stage[j] == test_value
        end
    end

    # Test final window spans all measurements
    n_windows = length(keys(prob_cl_ms.ms_windows))
    t_uni = PEtabTraining._get_unique_time_points(prob_original)
    @test t_uni[1] == prob_cl_ms.ms_windows[Symbol("stage$(n_windows)")][1][1]
    @test t_uni[end] == prob_cl_ms.ms_windows[Symbol("stage$(n_windows)")][1][end]

    # Test setting initial window values as Constant for each stage
    for i in 1:(n_windows - 1)
        prob_ms = prob_cl_ms.petab_problems[i]
        x_test = get_x(prob_ms)
        set_u0_ms_windows!(x_test, prob_cl_ms, i; init = MsInitConstant(5.0))

        u0_x_names = PEtabTraining._get_ms_u0_x_names(prob_ms)
        if window_u0_scale == :lin
            @test all(x_test[u0_x_names] .== 5.0)
        elseif window_u0_scale == :log10
            u0_x_names .= Symbol.("log10_" .* string.(u0_x_names))
            @test all(x_test[u0_x_names] .== log10(5.0))
        end
    end

    # Use initial values from first window
    x_original = get_x(prob_cl_ms.original)
    for i in 1:(n_windows - 1)
        prob_ms = prob_cl_ms.petab_problems[i]
        x_test = get_x(prob_ms)
        set_u0_ms_windows!(
            x_test, prob_cl_ms, i, x_original; init = MsInitFirst()
        )

        condition_ids = prob_ms.model_info.model.petab_tables[:conditions].conditionId
        for condition_id in condition_ids
            condition_id_original = PEtabTraining._get_condition_id_from_window(condition_id)

            u0_original = PEtab.get_u0(
                x_original, prob_original; retmap = false, condition = condition_id_original
            )
            u0_ms = PEtab.get_u0(x_test, prob_ms; retmap = false, condition = condition_id)
            @test all(.≈(u0_ms, u0_original; atol = 1.0e-8))
        end
    end

    # Test setting MS penalty
    for i in 1:(n_windows - 1)
        prob_ms = prob_cl_ms.petab_problems[i]
        x_test = get_x(prob_ms)
        if :net1 in keys(x_test)
            x_test.net1 .= 0.1
        end
        nllh1 = prob_ms.nllh(x_test)
        set_ms_window_penalty!(prob_cl_ms, 2.0)
        nllh2 = prob_ms.nllh(x_test)
        @test nllh1 != nllh2
        set_ms_window_penalty!(prob_cl_ms, 1.0)
    end

    # Check that an x-vector can be mapped between stages with util functions
    for i in 2:(n_windows - 1)
        x_from = deepcopy(get_x(prob_cl_ms.petab_problems[i - 1]))
        x_from .= 0.5
        x_to = PEtabTraining.map_x_stage(x_from, prob_cl_ms, i - 1, i)
        for (i, label) in pairs(ComponentArrays.labels(x_from))
            !(label in ComponentArrays.labels(x_to)) && continue
            ix = only(ComponentArrays.label2index(x_to, label))
            @test x_to[ix] == x_from[i]
        end
    end
    for i in n_windows:-1:2
        x_from = deepcopy(get_x(prob_cl_ms.petab_problems[i]))
        x_from .= 0.5
        x_to = PEtabTraining.map_x_stage(x_from, prob_cl_ms, i, i - 1)
        for (i, label) in pairs(ComponentArrays.labels(x_from))
            !(label in ComponentArrays.labels(x_to)) && continue
            ix = only(ComponentArrays.label2index(x_to, label))
            @test x_to[ix] == x_from[i]
        end
    end

    # Most easy by ensuring identical simulations between the stage and original problem.
    # See MS tests for comments
    for i in 1:(n_windows - 1)
        prob_ms = prob_cl_ms.petab_problems[i]
        x_original = get_x(prob_original)
        x_test = get_x(prob_ms)
        if :net1 in keys(x_test)
            x_original.net1 .= 0.001
            x_test.net1 .= 0.001
        end

        set_u0_ms_windows!(
            x_test, prob_cl_ms, i, x_original; init = MsInitSimulate()
        )

        ms_windows = prob_cl_ms.ms_windows[Symbol("stage$i")]
        prob_duplicated = _get_prob_duplicated(model_id, prob_original, ms_windows)

        nllh_cl_ms = prob_ms.nllh(x_test)
        nllh_cl_ms -= 0.5 * log(2π) * get_n_diff(prob_ms, prob_duplicated)
        nllh_duplicated = prob_duplicated.nllh(x_original)
        if model_id != "ude_model"
            @test nllh_cl_ms ≈ nllh_duplicated atol = 1.0e-3
        else
            @test nllh_cl_ms ≈ nllh_duplicated atol = 1.0e-2
        end
    end
    @test prob_cl_ms.petab_problems[end].nllh(x_original) == prob_original.nllh(x_original)
    return nothing
end

@testset "CL + MS" begin
    for n_windows in [2, 3, 5]
        test_cl_ms("Boehm_JProteomeRes2014", SplitTime(n_windows))
        test_cl_ms("mm_model_julia_defined", SplitTime(n_windows))
        test_cl_ms("ude_model", SplitTime(n_windows))
    end
    for n_windows in [2, 4]
        test_cl_ms("Fujita_SciSignal2010", SplitTime(n_windows))
    end

    # Test initial window parameters can be estimated on log10-scale
    test_cl_ms("mm_model_julia_defined", SplitTime(3); window_u0_scale = :log10)

    splits_test = [[15.0, 40.0, 100.0], [13.0, 25.0, 105.0]]
    for time_splits in splits_test
        test_cl_ms("Boehm_JProteomeRes2014", SplitTime(time_splits))
    end

    # Output regularization should be applied to each window for each stage
    prob_original = _get_petab_problem("ude_model"; include_regularization = true)
    cl_ms_prob = PEtabClMsProblem(
        prob_original, SplitTime(4); regularization_obs = "reg_o",
        regularization_specie = "nn_norm"
    )
    for prob in cl_ms_prob.petab_problems
        test_output_regularization_ms_prob(prob)
    end
end
