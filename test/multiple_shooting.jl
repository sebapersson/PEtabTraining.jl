using CSV, DataFrames, PEtab, PEtabTraining, Test

include(joinpath(@__DIR__, "helper.jl"))

function test_multiple_shooting(model_id, n_windows)
    path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
    prob_original = PEtabModel(path_yaml) |> PEtabODEProblem

    prob = PEtabMultipleShootingProblem(prob_original, SplitUniform(n_windows))
    a = 1
    mdf_original = prob_original.model_info.model.petab_tables[:measurements]
    mdf_ms = prob.petab_prob_ms.model_info.model.petab_tables[:measurements]

    # Check each window has correct start and end time-points. Note, the last measurement
    # point for a condition might appear in the middle of a window.
    unique_t = PEtabTraining._get_unique_timepoints(mdf_original)
    windows = PEtabTraining._makechunks(unique_t, n_windows; overlap = 1)
    cids = prob.petab_prob_ms.model_info.model.petab_tables[:conditions].conditionId |> unique
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
        m_cid_df_original = mdf_original[findall(x -> x ≥ tmin && x ≤ tmax, mdf_original.time), :]
        m_cid_df_original = m_cid_df_original[m_cid_df_original[!, :simulationConditionId] .== cid_original, :]
        next_cid = replace(cid, "__window$(window_index)_" => "__window$(window_index+1)_")
        if next_cid in cids
            @test nrow(m_cid_df_ms) == nrow(m_cid_df_original) + n_species
        else
            @test nrow(m_cid_df_ms) == nrow(m_cid_df_original)
        end
    end

    # Check the initial value parameters for each window can get properly set
    # Constant value
    xnames_u0 = PEtabTraining._get_ms_u0_xnames(prob)
    PEtabTraining.set_u0_windows!(prob, 4.0)
    @test all(prob.petab_prob_ms.xnominal[xnames_u0] .== 4.0)
    @test all(prob.petab_prob_ms.xnominal_transformed[xnames_u0] .== 4.0)
    # Initial values for the first window
    x_original = get_x(prob.original)
    PEtabTraining.set_u0_windows!(prob, x_original, :window1_u0)
    x_ms = get_x(prob.petab_prob_ms)
    for cid in cids
        cid_original = PEtabTraining._get_cid_from_window_id(cid)
        u0_original = PEtab.get_u0(x_original, prob.original; retmap = false, cid = cid_original)
        u0_ms = PEtab.get_u0(x_ms, prob.petab_prob_ms; retmap = false, cid = cid)
        @test u0_ms == u0_original
    end
    # Values from simulating initial values.
    # Easiest to check comparing the likelihood between original and ms problem, where for the
    # original problem, only difference is that duplicated points must be added. Note, due
    # to penalty entering via a likelihood, 0.5 * log(2π) must be accounted for.
    PEtabTraining.set_u0_windows!(prob, x_original, :window1_simulate)
    prob_duplicated = _get_prob_duplicated(prob.original, windows)
    nllh_ms = prob.petab_prob_ms.nllh(get_x(prob.petab_prob_ms))
    nllh_ms -= 0.5 * log(2π) * length(xnames_u0)
    nllh_duplicated = prob_duplicated.nllh(get_x(prob_duplicated))
    @test nllh_ms ≈ nllh_duplicated atol=1e-3
    return nothing
end

@testset "Multiple shooting" begin
    for n_windows in [2, 3, 5]
        test_multiple_shooting("Boehm_JProteomeRes2014", n_windows)
    end
    for n_windows in [2, 4]
        test_multiple_shooting("Fujita_SciSignal2010", n_windows)
    end
    # Test more windows is very heavy on RAM, due to huge number of parameters added
    test_multiple_shooting("Bachmann_MSB2011", 2)
end

# TODO: Error check pre-eq models not allowed
