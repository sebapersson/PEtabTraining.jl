struct PEtabClProblem{T <: Union{SplitTime, SplitData}}
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
    split_algorithm::T
    regularization_obs::Symbol
end
function PEtabClProblem(
        prob_original::PEtabODEProblem, split_alg::Union{SplitTime, SplitData};
        regularization_obs::Union{Nothing, String, Symbol} = nothing
    )::PEtabClProblem
    _check_regularization_obs(regularization_obs, prob_original)

    petab_tables_cl = _get_petab_tables_cl(split_alg, prob_original)

    _cl_adjust_ml_output_regularization!(petab_tables_cl, regularization_obs, prob_original)

    model_original = prob_original.model_info.model
    petab_problems = Vector{PEtabODEProblem}(undef, length(petab_tables_cl))
    for i in eachindex(petab_tables_cl)
        _filter_condition_table!(petab_tables_cl[i])

        if model_original.defined_in_julia == false
            model = PEtab._PEtabModel(
                model_original.paths, petab_tables_cl[i], false, false, true, false,
                model_original.petab_events, model_original.ml_models
            )
        else
            model = PEtab._PEtabModel(
                model_original.sys, petab_tables_cl[i], model_original.name,
                model_original.speciemap, model_original.parametermap,
                model_original.petab_events, model_original.ml_models, false;
                float_tspan = model_original.float_tspan
            )
        end
        petab_problems[i] = _PEtabODEProblem(model, prob_original)
    end

    regularization_obs = isnothing(regularization_obs) ? :none : regularization_obs
    return PEtabClProblem(petab_problems, prob_original, split_alg, regularization_obs)
end

function _get_petab_tables_cl(
        split_alg::Union{<:SplitTime, <:SplitData}, prob::PEtabODEProblem,
    )::Vector{PEtab.PEtabTables}
    measurements_df = _get_sorted_measurements_df(prob)

    chunks = _split_cl(split_alg, prob)
    out = Vector{PEtab.PEtabTables}(undef, length(chunks))
    for (i, chunk) in pairs(chunks)
        out[i] = copy(prob.model_info.model.petab_tables)
        if split_alg isa SplitTime
            out[i][:measurements] = measurements_df[measurements_df[!, :time] .≤ chunk, :]
        else
            out[i][:measurements] = measurements_df[1:chunk, :]
        end
    end
    return out
end

function _cl_adjust_ml_output_regularization!(
        petab_tables_cl, regularization_obs::Union{String, Symbol, Nothing},
        prob_original::PEtabODEProblem
    )::Nothing
    isnothing(regularization_obs) && return nothing

    regularization_obs = string(regularization_obs)
    measurements_original, conditions_original = PEtab._get_petab_tables(
        prob_original.model_info.model.petab_tables, [:measurements, :conditions]
    )

    for i in eachindex(petab_tables_cl)
        measurements_stage_i = petab_tables_cl[i][:measurements]
        for cid in conditions_original.conditionId
            measurements_stage_cid = filter(
                r -> r.simulationConditionId == cid, measurements_stage_i
            )
            if isempty(measurements_stage_cid)
                continue
            end

            # All conditions are not required to have output regularization
            measurements_cid_original = filter(
                r -> r.simulationConditionId == cid, measurements_original
            )
            if !(regularization_obs in measurements_cid_original.observableId)
                continue
            end

            if regularization_obs in measurements_stage_cid.observableId
                continue
            end

            row_reg = filter(
                r -> r.observableId == regularization_obs, measurements_cid_original
            )[1, :]
            row_reg.time = maximum(measurements_stage_cid.time)
            DataFrames.append!(
                measurements_stage_i, DataFrame(row_reg), promote = true, cols = :subset
            )
        end
    end
    return nothing
end
