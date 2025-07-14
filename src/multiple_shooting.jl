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
    model_ms = PEtab._PEtabModel(
        prob_original.model_info.model.paths, petab_tables_ms, false, false, true, false)

    @unpack (solver, solver_gradient, ss_solver, ss_solver_gradient, gradient_method, hessian_method, sensealg, reuse_sensitivities) = prob_original.probinfo
    petab_prob_ms = PEtabODEProblem(
        model_ms; odesolver = solver, odesolver_gradient = solver_gradient,
        ss_solver = ss_solver, ss_solver_gradient = ss_solver_gradient,
        gradient_method = gradient_method, hessian_method = hessian_method,
        sensealg = sensealg, reuse_sensitivities = reuse_sensitivities)
    # Set u0 values for each condition
    for cid in string.(petab_prob_ms.model_info.simulation_info.conditionids[:experiment])
        !occursin("__window", cid) && continue
        window_index = _get_index_from_window_id(cid)
        t_start = minimum(windows[window_index])
        petab_prob_ms.model_info.simulation_info.tstarts[Symbol(cid)] = t_start
    end
    return PEtabMultipleShootingProblem(split_algorithm, 1.0, petab_prob_ms, prob_original)
end

function _add_first_window!(
        petab_tables::Dict{Symbol, DataFrame}, condition_df_original::DataFrame,
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

    for cid in condition_df_original.conditionId
        _add_overlap_windows!(
            petab_tables, cid, condition_df_original, window, 1, specie_ids)
    end

    # Rename conditions to make it explicit which condition refers to the first window
    for cid in condition_df_original.conditionId
        i_row = findfirst(x -> x == cid, conditions_df.conditionId)
        _cid = _get_window_id(cid, 1, "", :condition)
        conditions_df[i_row, :conditionId] = _cid
    end

    _update_measurements!(petab_tables, condition_df_original, window, 1)
    return nothing
end

function _add_windows!(
        petab_tables::Dict{Symbol, DataFrame}, condition_df_original::DataFrame,
        speciemap::Vector, windows::Vector{<:Vector{<:Real}})::Nothing
    specie_ids = _get_specie_ids(speciemap)
    for i in 2:length(windows)
        for cid in condition_df_original.conditionId
            i == length(windows) && continue
            window = windows[i]
            _add_overlap_windows!(
                petab_tables, cid, condition_df_original, window, i, specie_ids)
        end
        _update_measurements!(petab_tables, condition_df_original, windows[i], i)
    end
    return nothing
end

function _add_window_penalty_parameter!(petab_tables::Dict{Symbol, DataFrame})::Nothing
    parameters_df = petab_tables[:parameters]
    df_ps = DataFrame(
        parameterId = "lambda_sqrt", parameterScale = "lin", lowerBound = 0.0,
        upperBound = Inf, nominalValue = 1.0, estimate = 0)
    append!(parameters_df, df_ps; promote = true, cols = :subset)
    return nothing
end

function _add_overlap_windows!(petab_tables::Dict{Symbol, DataFrame}, cid::String,
        condition_df_original::DataFrame, window::Vector{<:Real},
        window_index::Integer, specie_ids::Vector{String})::Nothing
    # If the condition does not have any measurements for time >maximum(window), no new
    # parameters need to be added, as the model should not be simulated past the window
    measurements_df = petab_tables[:measurements]
    i_rows = measurements_df[!, :simulationConditionId] .== cid
    if sum(measurements_df[i_rows, :time] .> maximum(window)) == 0
        return nothing
    end

    conditions_df = petab_tables[:conditions]
    parameters_df = petab_tables[:parameters]
    observable_df = petab_tables[:observables]

    # window_index + 1 needed, as setting initial conditions for the coming condition
    df_c = condition_df_original[condition_df_original.conditionId .== cid, :] |> deepcopy
    _cid = _get_window_id(cid, window_index + 1, "", :condition)
    df_c[1, :conditionId] = _cid
    for specie_id in specie_ids
        # Set initial value parameters
        pid = _get_window_id(cid, window_index, specie_id, :parameter)
        df_ps = DataFrame(parameterId = pid, parameterScale = "lin", lowerBound = 0.0,
            upperBound = Inf, nominalValue = 1e-3, estimate = 1)
        append!(parameters_df, df_ps; promote = true, cols = :subset)

        # Assign initial value to parameter value via conditions table
        if specie_id in names(df_c)
            df_c[i, specie_id] = pid
        else
            df_c[!, specie_id] .= pid
        end

        # Setup observables. Multiple window penalty is given by lambda_sqrt as the
        # parameter will be squared during the likelihood computations.
        obs_id = _get_window_id(cid, window_index, specie_id, :observable)
        df_obs = DataFrame(observableId = obs_id,
            observableFormula = "lambda_sqrt * ($(specie_id) - $(pid))",
            noiseFormula = "1.0", observableTransformation = "lin",
            noiseDistribution = "normal")
        append!(observable_df, df_obs, cols = :subset)

        # Add in the measurement table the penalty. Given the observable formula above,
        # the measurement value is set to zero. The different cid is needed, as the penalty
        # should be computed using the simulation from the preceding simulation.
        _cid_m = _get_window_id(cid, window_index, "", :condition)
        df_m = DataFrame(observableId = obs_id, simulationConditionId = _cid_m,
            measurement = 0.0, time = maximum(window))
        append!(measurements_df, df_m, cols = :subset)
    end
    append!(conditions_df, df_c)
    return nothing
end

function _update_measurements!(
        petab_tables::Dict{Symbol, DataFrame}, condition_df_original::DataFrame,
        window::Vector{<:Real}, window_index::Integer)::Nothing
    # Measurements can occur at the edge of a window, in this case observations must be
    # double counted in the measurement table and appear for the condition in the next
    # windows.
    measurements_df = petab_tables[:measurements]
    for cid in condition_df_original.conditionId
        im_cid = findall(x -> x == cid, measurements_df.simulationConditionId)
        im_time = findall(x -> x ≥ window[1] && x ≤ window[end], measurements_df.time)
        im_window = intersect(im_cid, im_time)
        _cid = _get_window_id(cid, window_index, "", :condition)
        measurements_df[im_window, :simulationConditionId] .= _cid

        # Not applicable to duplicate observations are on last window
        maximum(window) == maximum(measurements_df.time) && continue
        _im_overlap = findall(x -> x == maximum(window), measurements_df.time)
        im_overlap = intersect(_im_overlap, im_window)
        isnothing(im_overlap) && continue
        _cid = _get_window_id(cid, window_index + 1, "", :condition)
        overlap_df = measurements_df[im_overlap, :] |> deepcopy
        overlap_df[!, :simulationConditionId] .= _cid
        append!(measurements_df, overlap_df)
    end
    return nothing
end

function _get_window_id(
        cid::String, window_index::Integer, specie_id::String, which::Symbol)::String
    @assert which in [:condition, :parameter, :observable]
    if which == :condition
        return "__window$(window_index)_$(cid)__"
    elseif which == :parameter
        return "__window$(window_index)_$(cid)__$(specie_id)"
    elseif which == :observable
        return "__window$(window_index)_obs_$(cid)__$(specie_id)"
    end
end

function _get_cid_from_window_id(id::Union{String, Symbol})::String
    return string(match(r"__window\d+_([^_]+.*)__", "$id").captures[1])
end

function _get_index_from_window_id(id::Union{String, Symbol})::Int64
    return parse(Int64, match(r"window(\d+)", "$id").captures[1])
end

function _get_specie_id_from_window_id(id::Union{String, Symbol})::String
    return split("$id", "__")[end]
end
