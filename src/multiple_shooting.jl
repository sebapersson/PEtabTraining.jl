function PEtabMsProblem(
        prob_original::PEtabODEProblem, split_alg::SplitTime;
        regularization_obs::Union{Nothing, String, Symbol} = nothing,
        regularization_specie::Union{Nothing, String, Symbol} = nothing,
        window_u0_scale::Symbol = :lin, window_penalty_scale = :lin
    )::PEtabMsProblem
    if prob_original.model_info.simulation_info.has_pre_equilibration
        throw(ArgumentError("Multiple shooting is not supported for models with \
            pre-equilibration simulation conditions."))
    end

    @argcheck window_u0_scale in [:lin, :log, :log10]
    @argcheck window_penalty_scale in [:lin, :log, :log10]

    _check_regularization_specie(regularization_obs, regularization_specie)

    ms_windows = _split_ms(split_alg, prob_original)
    petab_ms_problem = _get_petab_ms_problem(
        prob_original, ms_windows, window_u0_scale, window_penalty_scale,
        _string(regularization_obs), _string(regularization_specie)
    )

    return PEtabMsProblem(
        petab_ms_problem, prob_original, split_alg, ms_windows, window_u0_scale,
        window_penalty_scale
    )
end

function _get_petab_ms_problem(
        prob_original::PEtabODEProblem, ms_windows::Vector{<:Vector{<:Real}},
        window_u0_scale::Symbol, window_penalty_scale::Symbol,
        regularization_obs::Union{Nothing, String},
        regularization_specie::Union{Nothing, String}
    )::PEtabODEProblem
    _check_regularization_obs(regularization_obs, prob_original)

    petab_tables_ms = deepcopy(prob_original.model_info.model.petab_tables)
    petab_tables_ms[:measurements] = DataFrame()

    # Values needed from original PEtabModel
    conditions_original, measurements_original = PEtab._get_petab_tables(
        prob_original.model_info.model.petab_tables, [:conditions, :measurements]
    )
    speciemap = prob_original.model_info.model.speciemap

    _add_first_window!(
        petab_tables_ms, measurements_original, conditions_original, ms_windows[1],
        speciemap, regularization_obs
    )
    _add_ms_windows!(
        petab_tables_ms, measurements_original, conditions_original, ms_windows,
        speciemap, window_u0_scale, window_penalty_scale, regularization_obs,
        regularization_specie
    )
    _add_window_penalty_parameter!(petab_tables_ms)

    _filter_condition_table!(petab_tables_ms)

    model_original = prob_original.model_info.model
    if model_original.defined_in_julia == false
        model_ms = PEtab._PEtabModel(
            prob_original.model_info.model.paths, petab_tables_ms, false, false, true,
            false, model_original.petab_events, prob_original.model_info.model.ml_models
        )
    else
        model_ms = PEtab._PEtabModel(
            model_original.sys, petab_tables_ms, model_original.name,
            model_original.speciemap, model_original.parametermap,
            model_original.petab_events, model_original.ml_models, false;
            float_tspan = model_original.float_tspan
        )
    end

    @unpack (
        solver, solver_gradient, ss_solver, ss_solver_gradient, gradient_method,
        hessian_method, sensealg, reuse_sensitivities,
    ) = prob_original.probinfo
    petab_ms_problem = PEtabODEProblem(
        model_ms; odesolver = solver, odesolver_gradient = solver_gradient,
        ss_solver = ss_solver, ss_solver_gradient = ss_solver_gradient,
        gradient_method = gradient_method, hessian_method = hessian_method,
        sensealg = sensealg, reuse_sensitivities = reuse_sensitivities
    )

    # In the PEtabODEProblem simulationInfo.tstarts must be altered, to ensure that each
    # simulations starts from the correct time-point to correctly handle potential events
    _adjust_tstarts!(petab_ms_problem, ms_windows)

    return petab_ms_problem
end

function _add_first_window!(
        petab_tables_ms::PEtab.PEtabTables, measurements_original::DataFrame,
        conditions_original::DataFrame, ms_window::Vector{<:Real},
        speciemap::Vector, regularization_obs
    )::Nothing
    specie_ids = _get_specie_ids(speciemap)
    conditions_df, parameters_df = PEtab._get_petab_tables(
        petab_tables_ms, [:conditions, :parameters]
    )

    # For the first window, initial values are already set in the PEtab problem, but
    # as all species will be assigned in the condition table, any specie not appearing
    # in the table must be moved into it
    for (i, specie_id) in pairs(specie_ids)
        specie_id in names(conditions_df) && continue

        # If initial value is an expression, NaN must be set to ensure the SBML
        # formula is used for the first window
        u0_value = string(speciemap[i].second)
        if !(PEtab.is_number(u0_value) || u0_value in parameters_df.parameterId)
            u0_value = "NaN"
        end
        conditions_df[!, specie_id] .= u0_value
    end

    # All original PEtab conditions must be re-named to explicitly map to a window
    for condition_id in conditions_original.conditionId
        row_idx = findfirst(x -> x == condition_id, conditions_df.conditionId)
        _condition_id = _get_window_id(condition_id, 1, "", :condition)
        conditions_df[row_idx, :conditionId] = _condition_id
    end

    _update_measurements!(
        petab_tables_ms, measurements_original, conditions_original, ms_window, 1,
        regularization_obs
    )
    return nothing
end

function _add_ms_windows!(
        petab_tables::PEtab.PEtabTables, measurements_original::DataFrame,
        conditions_original::DataFrame, ms_windows::Vector{<:Vector{<:Real}},
        speciemap::Vector, window_u0_scale::Symbol, window_penalty_scale::Symbol,
        regularization_obs, regularization_specie
    )::Nothing
    specie_ids = _get_specie_ids(speciemap)

    for i_window in 2:length(ms_windows)
        for condition_id in conditions_original.conditionId
            _add_overlap_ms_windows!(
                petab_tables, condition_id, measurements_original, conditions_original,
                ms_windows[i_window], i_window, specie_ids, window_u0_scale,
                window_penalty_scale, regularization_specie
            )
        end

        _update_measurements!(
            petab_tables, measurements_original, conditions_original, ms_windows[i_window],
            i_window, regularization_obs
        )
    end
    return nothing
end

function _add_overlap_ms_windows!(
        petab_tables_ms::PEtab.PEtabTables, condition_id::String,
        measurements_original::DataFrame, conditions_original::DataFrame,
        ms_window::Vector{<:Real}, i_window::Integer, specie_ids::Vector{String},
        window_u0_scale::Symbol, window_penalty_scale::Symbol, regularization_specie
    )::Nothing
    measurements_df, conditions_df, parameters_df, observable_df = PEtab._get_petab_tables(
        petab_tables_ms, [:measurements, :conditions, :parameters, :observables]
    )

    # If all time-points for the condition id are before the window, nothing to add
    measurements_condition_id = filter(
        row -> row.simulationConditionId == condition_id, measurements_original
    )
    if all(measurements_condition_id.time .< minimum(ms_window))
        return nothing
    end

    condition_id_window = _get_window_id(condition_id, i_window, "", :condition)
    condition_id_prev_window = _get_window_id(condition_id, i_window - 1, "", :condition)

    # Condition data-frame assigning initial values for the window
    condition_df = filter(row -> row.conditionId == condition_id, conditions_original) |>
        deepcopy
    condition_df[1, :conditionId] = condition_id_window
    for specie_id in specie_ids
        # Window initial value parameters
        parameter_id = _get_window_id(condition_id, i_window, specie_id, :parameter)

        if specie_id != regularization_specie
            lb = window_u0_scale == :lin ? -Inf : 1.0e-8
            df_ps = DataFrame(
                parameterId = parameter_id, parameterScale = "$(window_u0_scale)",
                lowerBound = lb, upperBound = 1.0e8, nominalValue = 1.0e-3, estimate = 1
            )
        else
            df_ps = DataFrame(
                parameterId = parameter_id, parameterScale = "lin", lowerBound = -Inf,
                upperBound = Inf, nominalValue = 0.0, estimate = 0
            )
        end
        DataFrames.append!(parameters_df, df_ps; promote = true, cols = :subset)

        # Assign the initial value for the window in condition table
        if specie_id in names(condition_df)
            condition_df[!, specie_id] .= string.(condition_df[!, specie_id])
            condition_df[1, specie_id] = parameter_id
        else
            condition_df[!, specie_id] .= parameter_id
        end

        # Output regularization species are not included in the penalty
        if specie_id == regularization_specie
            continue
        end

        # Observable. Multiple window penalty is given by lambda_sqrt as the
        # parameter will be squared during the likelihood computations.
        observable_id = _get_window_id(condition_id, i_window, specie_id, :observable)
        if window_penalty_scale == :lin
            obs_formula = "lambda_sqrt * ($(specie_id) - $(parameter_id))"
        else
            ft = "$window_penalty_scale"
            obs_formula = "lambda_sqrt * ($(ft)(abs($(specie_id))) - $(ft)($(parameter_id)))"
        end
        obs_df = DataFrame(
            observableId = observable_id, observableFormula = obs_formula,
            noiseFormula = "1.0", observableTransformation = "lin",
            noiseDistribution = "normal"
        )
        DataFrames.append!(observable_df, obs_df, cols = :subset)

        # Measurement for penalty
        measurements_row = DataFrame(
            observableId = observable_id, simulationConditionId = condition_id_prev_window,
            measurement = 0.0, time = minimum(ms_window)
        )
        if "simulationStartTime" in names(measurements_df)
            measurements_row[!, :simulationStartTime] .= 0.0
        end
        DataFrames.append!(measurements_df, measurements_row, cols = :subset)
    end
    DataFrames.append!(conditions_df, condition_df, promote = true)
    return nothing
end

function _add_window_penalty_parameter!(petab_tables_ms::PEtab.PEtabTables)::Nothing
    parameters_df = petab_tables_ms[:parameters]
    parameters_row = DataFrame(
        parameterId = "lambda_sqrt", parameterScale = "lin", lowerBound = 0.0,
        upperBound = 1.0e8, nominalValue = 1.0, estimate = 0
    )
    DataFrames.append!(parameters_df, parameters_row; promote = true, cols = :subset)
    return nothing
end

function _update_measurements!(
        petab_tables_ms::PEtab.PEtabTables, measurements_original::DataFrame,
        conditions_original::DataFrame, ms_window::Vector{<:Real},
        i_window::Integer, regularization_obs
    )::Nothing
    measurements_df = petab_tables_ms[:measurements]

    for condition_id in conditions_original.conditionId
        condition_id_window = _get_window_id(condition_id, i_window, "", :condition)

        idx_condition = findall(
            x -> x == condition_id, measurements_original.simulationConditionId
        )
        idx_time = findall(
            x -> x ≥ ms_window[1] && x ≤ ms_window[2], measurements_original.time
        )
        measurements_df_tmp = measurements_original[intersect(idx_condition, idx_time), :]

        _add_regularization_measurements!(
            measurements_df_tmp, regularization_obs, condition_id, measurements_original
        )

        measurements_df_tmp.simulationConditionId .= condition_id_window
        DataFrames.append!(
            measurements_df, measurements_df_tmp; cols = :subset, promote = true
        )
    end
    return nothing
end

function _get_window_id(
        condition_id::String, i_window::Integer, specie_id::String, which::Symbol
    )::String
    @assert which in [:condition, :parameter, :observable]

    if which == :condition
        return "_WINDOW$(i_window)_CONDITION_$(condition_id)"
    elseif which == :parameter
        return "_WINDOW$(i_window)_CONDITION_$(condition_id)_PARAMETER_$(specie_id)"
    else
        return "_WINDOW$(i_window)_CONDITION_$(condition_id)_OBSERVABLE_$(specie_id)"
    end
end

function _get_condition_id_from_window(s::AbstractString)::String
    m = match(r"_CONDITION_(.+?)(?=_(?:PARAMETER|OBSERVABLE)_|$)", s)
    @assert !isnothing(m) "Could not extract condition id from $s"
    return isnothing(m) ? nothing : m.captures[1]
end

function _get_index_from_window(id::Union{String, Symbol})::Int64
    return parse(Int64, match(r"_WINDOW(\d+)", "$id").captures[1])
end

function _get_specie_id_from_window(s::AbstractString)::String
    m = match(r"_(?:PARAMETER|OBSERVABLE)_(.+)$", s)
    @assert !isnothing(m) "Could not extract specie id from $s"
    return m.captures[1]
end

function _get_ms_u0_x_names(prob::PEtabODEProblem)::Vector{Symbol}
    names = Symbol[]
    for xname in prob.xnames
        !occursin(r"^_WINDOW\d+_CONDITION_", string(xname)) && continue
        push!(names, xname)
    end
    return names
end

function _add_regularization_measurements!(
        measurements_df_tmp::DataFrame, regularization_obs, condition_id::String,
        measurements_original::DataFrame
    )::Nothing
    isnothing(regularization_obs) && return nothing

    measurements_condition_original = filter(
        row -> row.simulationConditionId == condition_id, measurements_original
    )

    if (
            regularization_obs in measurements_condition_original.observableId &&
                regularization_obs in measurements_df_tmp.observableId
        )
        idx_row = findfirst(
            x -> x == regularization_obs, measurements_df_tmp.observableId
        )
        measurements_df_tmp.time[idx_row] = maximum(measurements_df_tmp.time)

    elseif regularization_obs in measurements_condition_original.observableId
        measurements_row = DataFrame(
            time = maximum(measurements_df_tmp.time), measurement = 0.0,
            observableId = regularization_obs, simulationConditionId = condition_id
        )

        if "simulationStartTime" in names(measurements_df_tmp)
            measurements_row[!, :simulationStartTime] .= 0.0
        end

        DataFrames.append!(
            measurements_df_tmp, measurements_row, promote = true, cols = :subset
        )
    end
    return nothing
end

function _adjust_tstarts!(
        petab_ms_problem::PEtabODEProblem, ms_windows::Vector{<:Vector{<:Real}}
    )::Nothing
    @unpack simulation_info = petab_ms_problem.model_info

    for condition_id in string.(simulation_info.conditionids[:experiment])
        !occursin("_WINDOW", condition_id) && continue
        i_window = _get_index_from_window(condition_id)
        t_start = ms_windows[i_window][1]
        simulation_info.tstarts[Symbol(condition_id)] = t_start
    end
    return nothing
end
