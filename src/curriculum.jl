struct PEtabCLProblem
    split_algorithm::Any
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
    regularization_obs::Symbol
end
function PEtabCLProblem(prob_original::PEtabODEProblem, split_algorithm; regularization_obs::Union{Nothing, String, Symbol} = nothing)::PEtabCLProblem
    _check_regularization_obs(regularization_obs, prob_original)

    model_original = prob_original.model_info.model
    petab_tables = _split(split_algorithm, prob_original, :curriculum)
    _cl_adjust_ml_output_regularization!(petab_tables, regularization_obs, prob_original)
    petab_problems = Vector{PEtabODEProblem}(undef, split_algorithm.nsplits)
    for i in eachindex(petab_tables)
        _filter_condition_table!(petab_tables[i])
        if prob_original.model_info.model.defined_in_julia == false
            model = PEtab._PEtabModel(model_original.paths, petab_tables[i], false,
                false, true, false, model_original.ml_models)
        else
            model = PEtab._PEtabModel(
                model_original.sys, petab_tables[i], model_original.name,
                model_original.speciemap, model_original.parametermap,
                model_original.callbacks, model_original.ml_models,
                false; float_tspan = model_original.float_tspan)
        end
        petab_problems[i] = _PEtabODEProblem(model, prob_original)
    end
    regularization_obs = isnothing(regularization_obs) ? :none : regularization_obs
    return PEtabCLProblem(split_algorithm, petab_problems, prob_original, regularization_obs)
end

function _split_curriculum(splits, mdf::DataFrame, prob::PEtabODEProblem, mode::Symbol)::Vector{PEtab.PEtabTables}
    @assert mode in [:datapoints, :time, :conditions]
    out = Vector{PEtab.PEtabTables}(undef, length(splits))
    if mode in [:datapoints, :time]
        for (i, split) in pairs(splits)
            out[i] = copy(prob.model_info.model.petab_tables)
            if mode == :time
                tmax = maximum(split)
                out[i][:measurements] = mdf[mdf[!, :time] .≤ tmax, :]
            else
                out[i][:measurements] = mdf[1:maximum(splits[i]), :]
            end
        end
    elseif mode == :conditions
        for (i, split) in pairs(splits)
            out[i] = copy(prob.model_info.model.petab_tables)
            irow = findall(x -> x in split, mdf.simulationConditionId)
            out[i][:measurements] = mdf[irow, :]
        end
    end
    return out
end

function _cl_adjust_ml_output_regularization!(petab_tables, regularization_obs::Union{String, Symbol, Nothing}, prob_original::PEtabODEProblem)::Nothing
    isnothing(regularization_obs) && return nothing
    regularization_obs = string(regularization_obs)

    measurements_original = prob_original.model_info.model.petab_tables[:measurements]
    conditions_original = prob_original.model_info.model.petab_tables[:conditions]
    for i in eachindex(petab_tables)
        measurements_stage = petab_tables[i][:measurements]
        for cid in conditions_original.conditionId
            measurements_stage_cid = filter(r -> r.simulationConditionId == cid, measurements_stage)
            if isempty(measurements_stage_cid)
                continue
            end

            # All conditions are not required to have output regularization
            measurements_cid_original = filter(r -> r.simulationConditionId == cid, measurements_original)
            if !(regularization_obs in measurements_cid_original.observableId)
                continue
            end

            if regularization_obs in measurements_stage_cid.observableId
                continue
            end

            row_reg = filter(r -> r.observableId == regularization_obs, measurements_cid_original)[1, :]
            row_reg.time = maximum(measurements_stage_cid.time)
            append!(measurements_stage, DataFrame(row_reg), promote = true, cols = :subset)
        end
    end
    return nothing
end
