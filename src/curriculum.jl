struct PEtabCLProblem
    split_algorithm::Any
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
end
function PEtabCLProblem(prob_original::PEtabODEProblem, split_algorithm)::PEtabCLProblem
    model_original = prob_original.model_info.model
    petab_tables = _split(split_algorithm, prob_original, :curriculum)
    petab_problems = Vector{PEtabODEProblem}(undef, split_algorithm.nsplits)
    for i in eachindex(petab_tables)
        _filter_condition_table!(petab_tables[i])
        if prob_original.model_info.model.defined_in_julia == false
            model = PEtab._PEtabModel(model_original.paths, petab_tables[i], false, false, true, false, model_original.ml_models)
        else
            model = PEtab._PEtabModel(model_original.sys, petab_tables[i], model_original.name, model_original.speciemap, model_original.parametermap, model_original.callbacks, model_original.ml_models, false; float_tspan = model_original.float_tspan)
        end
        petab_problems[i] = _PEtabODEProblem(model, prob_original)
    end
    return PEtabCLProblem(split_algorithm, petab_problems, prob_original)
end
function _split_curriculum(splits, mdf::DataFrame, prob::PEtabODEProblem,
        mode::Symbol)::Vector{PEtab.PEtabTables}
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
