include(joinpath(@__DIR__, "mm_model.jl"))
include(joinpath(@__DIR__, "ude_model.jl"))

function test_nllh(
        model_id, mdf::DataFrame, mdf_tmp::DataFrame, petab_prob, stage_problems_i)::Nothing
    if model_id == "mm_julia"
        petab_prob_ref = _get_mm_model(; measurements_df = mdf_tmp) |>
                         PEtabODEProblem
    elseif model_id == "ude"
        model = _get_lv_ude_model(; measurements_df = mdf_tmp)
        petab_prob_ref = PEtabODEProblem(model; odesolver = ODESolver(Rodas5P()))
    else
        CSV.write(petab_prob.model_info.model.paths[:measurements], mdf_tmp, delim = '\t')
        petab_prob_ref = _get_petab_problem(model_id)
        CSV.write(petab_prob.model_info.model.paths[:measurements], mdf, delim = '\t')
    end
    nllh_ref = petab_prob_ref.nllh(get_x(petab_prob_ref))
    nllh_test = stage_problems_i.nllh(get_x(petab_prob))
    if model_id != "ude"
        @test nllh_ref ≈ nllh_test atol=1e-8
    else
        @test nllh_ref ≈ nllh_test atol=1e-3
    end
    return nothing
end

function test_chunk_size(chunks)
    n_stages = length(chunks)
    for i in 1:(length(chunks) - 1)
        # Sometimes chunks cannot be divided uniformly, then the last is bigger
        if i == length(chunks) - 1
            @test abs(length(chunks[i]) - length(chunks[i + 1])) < n_stages
        else
            @test length(chunks[i]) == length(chunks[i + 1])
        end
    end
end

function _get_prob_duplicated(model_id, prob::PEtabODEProblem, ms_windows)
    model_original = prob.model_info.model
    mdf = prob.model_info.model.petab_tables[:measurements]
    mdf_duplicate = DataFrame()
    for window in ms_windows
        irow = findall(x -> x ≥ minimum(window) && x ≤ maximum(window), mdf.time)
        mdf_duplicate = vcat(mdf_duplicate, mdf[irow, :])
    end
    if model_id == "mm_julia"
        model = _get_mm_model(; measurements_df = mdf_duplicate)
    elseif model_id == "ude"
        model = _get_lv_ude_model(; measurements_df = mdf_duplicate)
    else
        tables_duplicate = deepcopy(prob.model_info.model.petab_tables)
        tables_duplicate[:measurements] = mdf_duplicate
        model = PEtab._PEtabModel(
            model_original.paths, tables_duplicate, false, false, true, false,
            model_original.petab_events, model_original.ml_models
        )
    end
    return PEtabODEProblem(model)
end

function _get_petab_problem(
        model_id::String; include_regularization::Bool = false
    )::PEtabODEProblem
    ode_solver = ODESolver(Rodas5P())
    if model_id == "mm_julia"
        model = _get_mm_model()
    elseif model_id == "ude"
        ode_solver = ODESolver(Rodas5P(), abstol=1e-12, reltol=1e-12)
        model = _get_lv_ude_model(; include_regularization = include_regularization)
    else
        path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
        model = PEtabModel(path_yaml)
    end
    return PEtabODEProblem(model; odesolver = ode_solver)
end

function get_n_diff(prob_ms::PEtabODEProblem, prob_duplicated::PEtabODEProblem)::Integer
    mdf_ms = prob_ms.model_info.model.petab_tables[:measurements]
    mdf_dup = prob_duplicated.model_info.model.petab_tables[:measurements]
    return nrow(mdf_ms) - nrow(mdf_dup)
end

function test_output_regularization_ms_prob(prob::PEtabODEProblem)::Nothing
    conditions = prob.model_info.model.petab_tables[:conditions]
    measurements = prob.model_info.model.petab_tables[:measurements]
    for condition_id in conditions.conditionId
        measurements_condition = filter(r -> r.simulationConditionId == condition_id, measurements)
        @test "reg_o" in measurements_condition.observableId
        i_rows = findall(x -> x == "reg_o", measurements_condition.observableId)
        @test length(i_rows) == 1
        @test measurements_condition.time[i_rows[1]] == maximum(measurements_condition.time)
    end
    # Sanity check nllh can be evaluated and ODE is correctly initialized
    @test !isinf(prob.nllh(get_x(prob)))
    for sol in values(prob.model_info.simulation_info.odesols)
        @test sol.prob.u0[end] == 0
    end
    return nothing
end
