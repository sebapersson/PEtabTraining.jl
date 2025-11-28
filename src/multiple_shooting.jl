mutable struct PEtabMultipleShootingProblem
    const split_algorithm::Any
    window_penalty::Float64
    const petab_prob_ms::PEtabODEProblem
    const original::PEtabODEProblem
end
function PEtabMultipleShootingProblem(prob_original::PEtabODEProblem,
        split_algorithm)::PEtabMultipleShootingProblem
    if prob_original.model_info.simulation_info.has_pre_equilibration
        throw(ArgumentError("Multiple shooting is not supported for models with \
            pre-equilibration simulation conditions."))
    end
    windows = _split(split_algorithm, prob_original, :multiple_shooting)
    petab_tables_ms = deepcopy(prob_original.model_info.model.petab_tables)

    # Values needed from original PEtabModel
    condition_df_original = prob_original.model_info.model.petab_tables[:conditions]
    speciemap = prob_original.model_info.model.speciemap

    # The first window is special, because each parameter already has an initial value,
    # which must be assigned in the condition table.
    PEtabTraining._add_first_window!(
        petab_tables_ms, condition_df_original, windows[1], speciemap)

    # Initial values for each window must be added as parameters, window penalty must be
    # added as observables + measurement points, and this must be done for each condition
    PEtabTraining._add_windows!(petab_tables_ms, condition_df_original, speciemap, windows)

    # Window penalty is added as a fixed parameters in the parameters table
    PEtabTraining._add_window_penalty_parameter!(petab_tables_ms)

    # For multiple-shooting, the entire setup can be captured in a singe PEtabModel and
    # subsequent PEtabODEProblem. In the PEtabODEProblem simulationInfo.tstarts must be
    # altered, to ensure that each simulations starts from the correct time-point to properly
    # handle any model events.
    _filter_condition_table!(petab_tables_ms)
    model_original = prob_original.model_info.model
    if model_original.defined_in_julia == false
        model_ms = PEtab._PEtabModel(prob_original.model_info.model.paths, petab_tables_ms, false, false, true, false, prob_original.model_info.model.ml_models)
    else
        model_ms = PEtab._PEtabModel(model_original.sys, petab_tables_ms, model_original.name, model_original.speciemap, model_original.parametermap, model_original.callbacks, model_original.ml_models, false; float_tspan = model_original.float_tspan)
    end

    @unpack (solver, solver_gradient, ss_solver, ss_solver_gradient, gradient_method,
        hessian_method, sensealg, reuse_sensitivities) = prob_original.probinfo
    petab_prob_ms = PEtabODEProblem(
        model_ms; odesolver = solver, odesolver_gradient = solver_gradient,
        ss_solver = ss_solver, ss_solver_gradient = ss_solver_gradient,
        gradient_method = gradient_method, hessian_method = hessian_method,
        sensealg = sensealg, reuse_sensitivities = reuse_sensitivities)
    # Set u0 values for each condition
    for cid in string.(petab_prob_ms.model_info.simulation_info.conditionids[:experiment])
        !occursin("__window", cid) && continue
        i_window = _get_index_from_window_id(cid)
        t_start = minimum(windows[i_window])
        petab_prob_ms.model_info.simulation_info.tstarts[Symbol(cid)] = t_start
    end
    return PEtabMultipleShootingProblem(split_algorithm, 1.0, petab_prob_ms, prob_original)
end

function _add_first_window!(
        petab_tables::PEtab.PEtabTables, condition_df_original::DataFrame,
        window::Vector{<:Real}, speciemap::Vector)::Nothing
    specie_ids = _get_specie_ids(speciemap)
    conditions_df = petab_tables[:conditions]
    parameters_df = petab_tables[:parameters]
    # Fix condition table for the first window (initial values already given)
    for i in eachindex(conditions_df.conditionId)
        for (j, specie_id) in pairs(specie_ids)
            specie_id in names(conditions_df) && continue
            # If initial value is an equation expression. In this case NaN must be set to
            # ensure the SBML formula is used for initial value computation
            u0_value = string(speciemap[j].second)
            if !(PEtab.is_number(u0_value) || u0_value in parameters_df.parameterId)
                u0_value = "NaN"
            end
            if specie_id in names(conditions_df)
                conditions_df[i, specie_id] = u0_value
            else
                conditions_df[!, specie_id] .= u0_value
            end
        end
    end

    # All windows must be re-named to explicitly map to a window
    for cid in condition_df_original.conditionId
        i_row = findfirst(x -> x == cid, conditions_df.conditionId)
        _cid = _get_window_id(cid, 1, "", :condition)
        conditions_df[i_row, :conditionId] = _cid
    end

    _update_measurements!(petab_tables, condition_df_original, window, 1)
    return nothing
end

function _add_windows!(petab_tables::PEtab.PEtabTables, condition_df_original::DataFrame, speciemap::Vector, windows::Vector{<:Vector{<:Real}})::Nothing
    specie_ids = _get_specie_ids(speciemap)
    for i_window in 2:length(windows)
        for cid in condition_df_original.conditionId
            _add_overlap_windows!(petab_tables, cid, condition_df_original, windows[i_window], i_window, specie_ids)
        end
        _update_measurements!(petab_tables, condition_df_original, windows[i_window], i_window)
    end
    return nothing
end

function _add_window_penalty_parameter!(petab_tables::PEtab.PEtabTables)::Nothing
    parameters_df = petab_tables[:parameters]
    df_ps = DataFrame(
        parameterId = "lambda_sqrt", parameterScale = "lin", lowerBound = 0.0,
        upperBound = Inf, nominalValue = 1.0, estimate = 0)
    append!(parameters_df, df_ps; promote = true, cols = :subset)
    return nothing
end

function _add_overlap_windows!(petab_tables::PEtab.PEtabTables, cid::String, condition_df_original::DataFrame, window::Vector{<:Real}, i_window::Integer, specie_ids::Vector{String})::Nothing
    measurements_df = petab_tables[:measurements]
    conditions_df = petab_tables[:conditions]
    parameters_df = petab_tables[:parameters]
    observable_df = petab_tables[:observables]

    # If all time-points for the condition id are before the window, no need to add
    measurements_cid_df = filter(row -> row.simulationConditionId == cid, measurements_df)
    if all(measurements_cid_df.time .< minimum(window))
        return nothing
    end

    df_cid = condition_df_original[condition_df_original.conditionId .== cid, :] |> deepcopy
    cid_window = _get_window_id(cid, i_window, "", :condition)
    cid_prev_window = _get_window_id(cid, i_window - 1, "", :condition)
    df_cid[1, :conditionId] = cid_window
    for specie_id in specie_ids
        # Create initial PEtab parameters for the window
        pid = _get_window_id(cid, i_window, specie_id, :parameter)
        df_ps = DataFrame(parameterId = pid, parameterScale = "lin", lowerBound = 0.0,
            upperBound = Inf, nominalValue = 1e-3, estimate = 1)
        append!(parameters_df, df_ps; promote = true, cols = :subset)

        # Assign initial value to parameter value via conditions table
        if specie_id in names(df_cid)
            df_cid[!, specie_id] = string.(df_cid[!, specie_id])
            df_cid[1, specie_id] = pid
        else
            df_cid[!, specie_id] .= pid
        end

        # Observable. Multiple window penalty is given by lambda_sqrt as the
        # parameter will be squared during the likelihood computations.
        obs_id = _get_window_id(cid, i_window, specie_id, :observable)
        df_obs = DataFrame(observableId = obs_id,
            observableFormula = "lambda_sqrt * ($(specie_id) - $(pid))",
            noiseFormula = "1.0", observableTransformation = "lin",
            noiseDistribution = "normal")
        append!(observable_df, df_obs, cols = :subset)

        # NOTE the cid (value computed from simulating previous window)
        df_m = DataFrame(observableId = obs_id, simulationConditionId = cid_prev_window,
            measurement = 0.0, time = minimum(window))
        append!(measurements_df, df_m, cols = :subset)
    end
    append!(conditions_df, df_cid, promote = true)
    return nothing
end

function _update_measurements!(
        petab_tables::PEtab.PEtabTables, condition_df_original::DataFrame,
        window::Vector{<:Real}, i_window::Integer)::Nothing
    # Measurements can occur at the edge of a window, in this case observations must be
    # double counted in the measurement table and appear for the condition in the next
    # windows.
    measurements_df = petab_tables[:measurements]
    for cid in condition_df_original.conditionId
        im_cid = findall(x -> x == cid, measurements_df.simulationConditionId)
        im_time = findall(x -> x ≥ window[1] && x ≤ window[end], measurements_df.time)
        im_window = intersect(im_cid, im_time)
        _cid = _get_window_id(cid, i_window, "", :condition)
        measurements_df[im_window, :simulationConditionId] .= _cid

        # Not applicable to duplicate observations are on last window
        maximum(window) == maximum(measurements_df.time) && continue
        _im_overlap = findall(x -> x == maximum(window), measurements_df.time)
        im_overlap = intersect(_im_overlap, im_window)
        isnothing(im_overlap) && continue
        _cid = _get_window_id(cid, i_window + 1, "", :condition)
        overlap_df = measurements_df[im_overlap, :] |> deepcopy
        overlap_df[!, :simulationConditionId] .= _cid
        append!(measurements_df, overlap_df)
    end
    return nothing
end

function _get_window_id(
        cid::String, i_window::Integer, specie_id::String, which::Symbol)::String
    @assert which in [:condition, :parameter, :observable]
    if which == :condition
        return "___window$(i_window)___$(cid)___"
    elseif which == :parameter
        return "___window$(i_window)___$(cid)___$(specie_id)"
    elseif which == :observable
        return "___window$(i_window)___obs___$(cid)___$(specie_id)"
    end
end

function _get_cid_from_window_id(id::Union{String, Symbol})::String
    s = string(id)
    m = match(r"^___window\d+___(.*)___.*$", s)
    return m.captures[1]
end

function _get_index_from_window_id(id::Union{String, Symbol})::Int64
    return parse(Int64, match(r"window(\d+)", "$id").captures[1])
end

function _get_specie_id_from_window_id(id::Union{String, Symbol})::String
    s = string(id)
    m = match(r"^___window\d+___(.*)___(.*)$", s)
    return m.captures[2]
end
