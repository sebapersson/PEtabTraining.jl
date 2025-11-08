using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "helper.jl"))

function test_multiple_shooting(model_id, n_windows::Integer)
    split_algorithm = SplitUniform(n_windows)
    _test_multiple_shooting(model_id, split_algorithm)
    return nothing
end
function test_multiple_shooting(model_id, windows::Vector)
    split_algorithm = SplitCustom(windows; mode = :time)
    _test_multiple_shooting(model_id, split_algorithm)
    return nothing
end

function _test_multiple_shooting(model_id, split_algorithm)
    prob_original = _get_petab_problem(model_id)

    prob = PEtabMultipleShootingProblem(prob_original, split_algorithm)

    mdf_original = prob_original.model_info.model.petab_tables[:measurements]
    mdf_ms = prob.petab_prob_ms.model_info.model.petab_tables[:measurements]

    # Check each window has correct start and end time-points. Note, the last measurement
    # point for a condition might appear in the middle of a window.
    unique_t = PEtabTraining._get_unique_timepoints(mdf_original)
    if split_algorithm isa SplitUniform
        windows = PEtabTraining._makechunks(unique_t, split_algorithm.nsplits; overlap = 1)
    else
        windows = PEtabTraining._splits_to_windows(split_algorithm.splits)
    end
    cids = prob.petab_prob_ms.model_info.model.petab_tables[:conditions].conditionId |>
           unique
    for cid in cids
        window_index = PEtabTraining._get_index_from_window_id(cid)
        t_min = minimum(windows[window_index])
        t_max = mdf_ms[mdf_ms[!, :simulationConditionId] .== cid, :].time |> maximum
        @test t_min == prob.petab_prob_ms.model_info.simulation_info.tstarts[Symbol(cid)]
        @test t_max == prob.petab_prob_ms.model_info.simulation_info.tmaxs[Symbol(cid)]
    end

    # Check each window has correct number of measurement points. As each window starts and
    # ends at data-points, these data-points should be double-counted (hence the ≥), and
    # each penalty gets a data-point, unless a condition does not have data-points for
    # the next window
    n_species = length(prob_original.model_info.model.speciemap)
    for cid in cids
        window_index = PEtabTraining._get_index_from_window_id(cid)
        tmin, tmax = minimum(windows[window_index]), maximum(windows[window_index])
        cid_original = PEtabTraining._get_cid_from_window_id(cid)
        m_cid_df_ms = mdf_ms[mdf_ms[!, :simulationConditionId] .== cid, :]
        m_cid_df_original = mdf_original[
        findall(x -> x ≥ tmin && x ≤ tmax, mdf_original.time), :]
        m_cid_df_original = m_cid_df_original[
        m_cid_df_original[
            !, :simulationConditionId] .== cid_original, :]
        next_cid = replace(cid, "__window$(window_index)_" => "__window$(window_index+1)_")
        if next_cid in cids
            @test nrow(m_cid_df_ms) == nrow(m_cid_df_original) + n_species
        else
            @test nrow(m_cid_df_ms) == nrow(m_cid_df_original)
        end
    end

    # Check the initial value parameters for each window can get properly set
    # - Constant value
    xnames_u0 = PEtabTraining._get_ms_u0_xnames(prob)
    PEtabTraining.set_u0_windows!(prob, 4.0)
    @test all(prob.petab_prob_ms.xnominal[xnames_u0] .== 4.0)
    @test all(prob.petab_prob_ms.xnominal_transformed[xnames_u0] .== 4.0)
    # - Initial values for the first window
    x_original = get_x(prob.original)
    PEtabTraining.set_u0_windows!(prob, x_original, :window1_u0)
    x_ms = get_x(prob.petab_prob_ms)
    for cid in cids
        cid_original = PEtabTraining._get_cid_from_window_id(cid)
        u0_original = PEtab.get_u0(
            x_original, prob.original; retmap = false, cid = cid_original)
        u0_ms = PEtab.get_u0(x_ms, prob.petab_prob_ms; retmap = false, cid = cid)
        @test u0_ms == u0_original
    end

    # Test effect changing penalty parameter (likelihood should change)
    # - Constant value
    nllh1 = prob.petab_prob_ms.nllh(x_ms)
    PEtabTraining.set_window_penalty!(prob, 2.0)
    nllh2 = prob.petab_prob_ms.nllh(x_ms)
    @test nllh1 != nllh2
    @test prob.petab_prob_ms.model_info.petab_parameters.nominal_value[end] ≈ √2
    # - Simulating from initial values
    # Easiest to check comparing the likelihood between original and ms problem, where for
    # the original problem, only difference is that duplicated points must be added. Note,
    # due to penalty entering via a likelihood, 0.5 * log(2π) must be accounted for. Not yet
    # applicable for provided as ODEProblem
    if prob.original.model_info.model.sys isa ODEProblem
        return nothing
    end
    PEtabTraining.set_u0_windows!(prob, x_original, :window1_simulate)
    prob_duplicated = _get_prob_duplicated(model_id, prob.original, windows)
    nllh_ms = prob.petab_prob_ms.nllh(get_x(prob.petab_prob_ms))
    nllh_ms -= 0.5 * log(2π) * length(xnames_u0)
    nllh_duplicated = prob_duplicated.nllh(get_x(prob_duplicated))
    @test nllh_ms≈nllh_duplicated atol=1e-3
    return nothing
end

function test_reference()
    prob_original = _get_petab_problem("mm_julia")
    prob_original.probinfo.solver.abstol = 1e-12
    prob_original.probinfo.solver.reltol = 1e-12
    prob_original.probinfo.solver.maxiters = Int(1e6)
    prob_ms = PEtabMultipleShootingProblem(prob_original, SplitCustom([3.25, 5.25, 10.0]; mode = :time))
    PEtabTraining.set_u0_windows!(prob_ms, 10.0)

    measurements_df = prob_ms.petab_prob_ms.model_info.model.petab_tables[:measurements]
    ix_not_window = findall(.!startswith.(measurements_df.observableId, "___window"))
    x_original = get_x(prob_original)
    x_ms = get_x(prob_ms.petab_prob_ms)

    # Test a naive computation of MS loss works
    residuals = prob_ms.petab_prob_ms.residuals(x_ms)
    nllh_ms = prob_ms.petab_prob_ms.nllh(x_ms)
    nllh_ms_naive = sum(@. 0.5 * log(2π) + 0.5 * residuals^2)
    @test nllh_ms ≈ nllh_ms_naive atol=1e-8

    # Check original loss can be retrieved
    PEtabTraining.set_u0_windows!(prob_ms, x_original, :window1_simulate)
    x_ms = get_x(prob_ms.petab_prob_ms)
    nllh_original = prob_ms.original.nllh(x_original)
    residuals = prob_ms.petab_prob_ms.residuals(x_ms)
    nllh_original_naive = sum(@. 0.5 * log(2π) + 0.5 * residuals[ix_not_window]^2)
    @test nllh_original ≈ nllh_original_naive atol=1e-6
    return nothing
end

@testset "Multiple shooting" begin
    for n_windows in [2, 3, 5]
        test_multiple_shooting("Boehm_JProteomeRes2014", n_windows)
        test_multiple_shooting("mm_julia", n_windows)
        test_multiple_shooting("ude", n_windows)
    end
    splits_test = [[15.0, 40.0, 100.0, 240.0], [13.0, 25.0, 105.0, 250.0]]
    for split in splits_test
        test_multiple_shooting("Boehm_JProteomeRes2014", split)
    end
    test_multiple_shooting("mm_julia", [5.0, 10.0])

    for n_windows in [2, 4]
        test_multiple_shooting("Fujita_SciSignal2010", n_windows)
    end
    # Test more windows is very heavy on RAM, due to huge number of parameters added
    test_multiple_shooting("Bachmann_MSB2011", 2)
    # Reference with manually computed ms-loss
    test_reference()
end
