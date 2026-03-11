using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "common.jl"))

function test_multiple_shooting(model_id, split_alg::SplitTime; window_u0_scale = :lin)
    @assert window_u0_scale in [:lin, :log10]

    prob_original = _get_petab_problem(model_id)

    prob_ms = PEtabMsProblem(prob_original, split_alg; window_u0_scale = window_u0_scale)

    measurements_original = prob_original.model_info.model.petab_tables[:measurements]
    measurements_ms = prob_ms.petab_ms_problem.model_info.model.petab_tables[:measurements]

    # Check each window has correct start and end time-points. Note, the last measurement
    # point for a condition might appear in the middle of a window.
    unique_t = PEtabTraining._get_unique_time_points(prob_original)
    if split_alg.spec isa Integer
        ms_windows = PEtabTraining._makechunks(unique_t, split_alg.spec; overlap = 1)
    else
        vals = [unique_t[1], split_alg.spec..., unique_t[end]]
        ms_windows = [[vals[i], vals[i + 1]] for i in 1:(length(vals) - 1)]
    end

    condition_ids = prob_ms.petab_ms_problem.model_info.model.petab_tables[:conditions].conditionId |>
        unique

    # Check t0 and t_end values are correct for each condition window
    for condition_id in condition_ids
        i_window = PEtabTraining._get_index_from_window(condition_id)
        t_min = minimum(ms_windows[i_window])
        idx_t = measurements_ms[!, :simulationConditionId] .== condition_id
        t_max = maximum(measurements_ms.time[idx_t])

        simulation_info = prob_ms.petab_ms_problem.model_info.simulation_info
        @test t_min == simulation_info.tstarts[Symbol(condition_id)]
        @test t_max == simulation_info.tmaxs[Symbol(condition_id)]
    end

    # Check each window has correct number of measurement points. As each window starts and
    # ends at data-points, these data-points should be double-counted (hence the ≥), and
    # each penalty gets a data-point, unless a condition does not have data-points for
    # the next window
    n_species = length(prob_original.model_info.model.speciemap)
    for condition_id in condition_ids
        i_window = PEtabTraining._get_index_from_window(condition_id)
        t_min = first(ms_windows[i_window])
        t_max = last(ms_windows[i_window])

        condition_id_original = PEtabTraining._get_condition_id_from_window(condition_id)
        measurements_original_condition = filter(
            row -> row.simulationConditionId == condition_id_original && row.time ≥ t_min &&
                row.time ≤ t_max, measurements_original
        )
        measurements_ms_condition = filter(
            row -> row.simulationConditionId == condition_id && row.time ≥ t_min &&
                row.time ≤ t_max, measurements_ms
        )

        next_condition_id = replace(
            condition_id, "_WINDOW$(i_window)_" => "_WINDOW$(i_window + 1)_"
        )
        if next_condition_id in condition_ids
            @test nrow(measurements_ms_condition) == (
                nrow(measurements_original_condition) + n_species
            )
        else
            @test nrow(measurements_ms_condition) == nrow(measurements_original_condition)
        end
    end

    # Check the initial value parameters for each window can get properly set
    # - Constant value
    u0_x_names = PEtabTraining._get_ms_u0_x_names(prob_ms.petab_ms_problem)
    x_test = get_x(prob_ms.petab_ms_problem)
    set_u0_ms_windows!(x_test, prob_ms; init = MsInitConstant(4.0))
    if window_u0_scale == :lin
        @test all(x_test[u0_x_names] .== 4.0)
    elseif window_u0_scale == :log10
        u0_x_names .= Symbol.("log10_" .* string.(u0_x_names))
        @test all(x_test[u0_x_names] .== log10(4.0))
    end

    # - Initial values for the first window
    x_original = get_x(prob_ms.original)
    x_test = get_x(prob_ms.petab_ms_problem)
    set_u0_ms_windows!(x_test, prob_ms, x_original; init = MsInitFirst())
    for condition_id in condition_ids
        condition_id_original = PEtabTraining._get_condition_id_from_window(condition_id)
        u0_original = PEtab.get_u0(
            x_original, prob_ms.original; retmap = false, condition = condition_id_original
        )
        u0_ms = PEtab.get_u0(
            x_test, prob_ms.petab_ms_problem; retmap = false, condition = condition_id
        )
        @test all(.≈(u0_ms, u0_original; atol = 1.0e-8))
    end

    # Test effect changing penalty parameter (likelihood should change)
    x_ms = get_x(prob_ms.petab_ms_problem)
    if :net1 in keys(x_ms)
        x_ms.net1 .= 0.01
    end
    nllh1 = prob_ms.petab_ms_problem.nllh(x_ms)
    set_ms_window_penalty!(prob_ms, 2.0)
    nllh2 = prob_ms.petab_ms_problem.nllh(x_ms)
    @test nllh1 != nllh2
    @test prob_ms.petab_ms_problem.model_info.petab_parameters.nominal_value[end] ≈ √2

    x_original = get_x(prob_ms.original)
    x_test = get_x(prob_ms.petab_ms_problem)
    if :net1 in keys(x_test)
        x_original.net1 .= 0.001
        x_test.net1 .= 0.001
    end
    set_u0_ms_windows!(x_test, prob_ms, x_original; init = MsInitSimulate())
    prob_duplicated = _get_prob_duplicated(model_id, prob_ms.original, ms_windows)

    nllh_ms = prob_ms.petab_ms_problem.nllh(x_test)
    nllh_ms -= 0.5 * log(2π) * length(u0_x_names)
    nllh_duplicated = prob_duplicated.nllh(x_original)
    # For UDE magnitude of likelihood is on order of 1e6, so numerics play a large role
    if model_id != "ude_model"
        @test nllh_ms ≈ nllh_duplicated atol = 1.0e-3
    else
        @test nllh_ms ≈ nllh_duplicated atol = 1.0e-2
    end
    return nothing
end

function test_reference()
    prob_original = _get_petab_problem("mm_model_julia_defined")
    prob_original.probinfo.solver.abstol = 1.0e-12
    prob_original.probinfo.solver.reltol = 1.0e-12
    prob_original.probinfo.solver.maxiters = Int(1.0e6)
    prob_ms = PEtabMsProblem(prob_original, SplitTime([3.25, 5.25]))

    x_test = get_x(prob_ms.petab_ms_problem)
    set_u0_ms_windows!(x_test, prob_ms; init = MsInitConstant(10.0))

    # Test a naive computation of MS loss works
    residuals = prob_ms.petab_ms_problem.residuals(x_test)
    nllh_ms = prob_ms.petab_ms_problem.nllh(x_test)
    nllh_ms_naive = sum(@. 0.5 * log(2π) + 0.5 * residuals^2)
    @test nllh_ms ≈ nllh_ms_naive atol = 1.0e-8
    return nothing
end

@testset "Multiple shooting" begin
    for n_windows in [2, 3, 5]
        test_multiple_shooting("Boehm_JProteomeRes2014", SplitTime(n_windows))
        test_multiple_shooting("mm_model_julia_defined", SplitTime(n_windows))
        test_multiple_shooting("ude_model", SplitTime(n_windows))
    end
    # This model has fewer data-points than the above
    for n_windows in [2, 4]
        test_multiple_shooting("Fujita_SciSignal2010", SplitTime(n_windows))
    end

    # Test that estimating window-parameters on log10-scale is feasible
    test_multiple_shooting("mm_model_julia_defined", SplitTime(3); window_u0_scale = :log10)

    splits_test = [[15.0, 40.0, 100.0], [13.0, 25.0, 105.0]]
    for time_splits in splits_test
        test_multiple_shooting("Boehm_JProteomeRes2014", SplitTime(time_splits))
    end
    test_multiple_shooting("mm_model_julia_defined", SplitTime([5.0, 9.0]))

    # Reference with manually computed ms-loss
    test_reference()

    # Output regularization should be applied to each window
    prob_original = _get_petab_problem("ude_model"; include_regularization = true)
    prob_ms = PEtabMsProblem(
        prob_original, SplitTime(4); regularization_obs = "reg_o",
        regularization_specie = "nn_norm"
    )
    test_output_regularization_ms_prob(prob_ms.petab_ms_problem)

    # Test more windows is very heavy on RAM, due to huge number of parameters added to
    # the model, but it is a good test-case
    test_multiple_shooting("Bachmann_MSB2011", SplitTime(2))
end
