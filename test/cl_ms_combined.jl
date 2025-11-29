using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "helper.jl"))

function test_ml_cl(model_id, n_splits)
    prob_original = _get_petab_problem(model_id)
    prob_cl_ms = PEtabCLMSProblem(prob_original, SplitUniform(n_splits))

    # Test that windows are created correctly
    mdf = prob_original.model_info.model.petab_tables[:measurements]
    t_unique = PEtabTraining._get_unique_timepoints(mdf)
    for i in length(prob_cl_ms.windows):-1:2
        windows_stage = prob_cl_ms.windows[Symbol("stage$i")]
        windows_stage_prev = prob_cl_ms.windows[Symbol("stage$(i-1)")]
        for j in eachindex(windows_stage)
            @test windows_stage[j] == unique(reduce(vcat, windows_stage_prev[j:(j+1)]))
        end
    end
    n_windows = length(keys(prob_cl_ms.windows))
    @test t_unique == prob_cl_ms.windows[Symbol("stage$(n_windows)")][1]

    # Test setting initial window values
    # Constant
    PEtabTraining.set_u0_windows!(prob_cl_ms, :constant, 4.0)
    for i in 1:(n_windows-1)
        _prob = prob_cl_ms.petab_problems[i]
        xnames_u0 = PEtabTraining._get_ms_u0_xnames(_prob)
        @test all(_prob.xnominal[xnames_u0] .== 4.0)
        @test all(_prob.xnominal_transformed[xnames_u0] .== 4.0)
    end
    # Initial values for the first window
    x_original = get_x(prob_cl_ms.original)
    PEtabTraining.set_u0_windows!(prob_cl_ms, :window1_u0, x_original)
    for i in 1:(n_windows-1)
        _prob = prob_cl_ms.petab_problems[i]
        x_ms = get_x(_prob)
        cids = _prob.model_info.model.petab_tables[:conditions].conditionId |> unique
        for cid in cids
            cid_original = PEtabTraining._get_cid_from_window_id(cid)
            u0_original = PEtab.get_u0(x_original, prob_original; retmap = false, cid = cid_original)
            u0_ms = PEtab.get_u0(x_ms, _prob; retmap = false, cid = cid)
            @test u0_ms == u0_original
        end
    end

    # Test setting MS penalty
    for i in 1:(n_windows-1)
        _prob = prob_cl_ms.petab_problems[i]
        x_ms = get_x(_prob)
        nllh1 = _prob.nllh(x_ms)
        PEtabTraining.set_window_penalty!(prob_cl_ms, 2.0)
        nllh2 = _prob.nllh(x_ms)
        @test nllh1 != nllh2
        PEtabTraining.set_window_penalty!(prob_cl_ms, 1.0)
    end

    # Check that an x-vector can be mapped between stages with util functions
    for i in 2:(n_windows-1)
        x_from = deepcopy(get_x(prob_cl_ms.petab_problems[i-1]))
        x_from .= 0.5
        x_to = PEtabTraining.map_x_stage(x_from, prob_cl_ms, i-1, i)
        for (i, label) in pairs(ComponentArrays.labels(x_from))
            !(label in ComponentArrays.labels(x_to)) && continue
            ix = only(ComponentArrays.label2index(x_to, label))
            @test x_to[ix] == x_from[i]
        end
    end
    for i in n_windows:-1:2
        x_from = deepcopy(get_x(prob_cl_ms.petab_problems[i]))
        x_from .= 0.5
        x_to = PEtabTraining.map_x_stage(x_from, prob_cl_ms, i, i-1)
        for (i, label) in pairs(ComponentArrays.labels(x_from))
            !(label in ComponentArrays.labels(x_to)) && continue
            ix = only(ComponentArrays.label2index(x_to, label))
            @test x_to[ix] == x_from[i]
        end
    end

    # Check correctness. Most easy by ensuing identical simulations between the stage and
    # original problem. See MS tests for comments
    if prob_original.model_info.model.sys isa ODEProblem
        return nothing
    end
    for i in 1:(n_windows-1)
        x_original = get_x(prob_original)
        _prob = prob_cl_ms.petab_problems[i]
        windows = prob_cl_ms.windows[Symbol("stage$i")]
        PEtabTraining.set_u0_windows!(prob_cl_ms, :window1_simulate, x_original)
        prob_duplicated = _get_prob_duplicated(model_id, prob_original, windows)
        nllh_cl_ms = _prob.nllh(get_x(_prob))
        nllh_cl_ms -= 0.5 * log(2π) * get_n_diff(_prob, prob_duplicated)
        nllh_duplicated = prob_duplicated.nllh(x_original)
        @test nllh_cl_ms≈nllh_duplicated atol=1e-3
    end
    @test prob_cl_ms.petab_problems[end].nllh(x_original) == prob_original.nllh(x_original)
    return nothing
end

@testset "CL + MS" begin
    for n in [2, 3, 5]
        test_ml_cl("Boehm_JProteomeRes2014", n)
        test_ml_cl("mm_julia", n)
        test_ml_cl("ude", n)
    end
    for n in [2, 4]
        test_ml_cl("Fujita_SciSignal2010", n)
    end
end
