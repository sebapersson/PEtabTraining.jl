using CSV, DataFrames, PEtab, PEtabTraining, Test

# This is something that is implemented at the PEtab-level, adding a huge new bunch
# of observables, parameters, and measurements in the measurements table and create
# entirely new simulations conditions. Further, must allow PEtab.jl to set its own
# start time, otherwise, events are not going to work.
# Naturally, a question becomes, how to deal with initial conditions for the base
# condition where the initial value might be an expression. The only natural way
# to deal with this is via a NaN check, where if the parameter setting the initial
# value is NaN, resorts the the SBML formula. This might cause gradient problems,
# not really as we never need to compile a trace. Similar, for Julia interface
# I do not really see any problem. The drawback, this will be a major pain
# to implement.

model_id = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "published_models", model_id, "$(model_id).yaml")
petab_prob = PEtabModel(path_yaml) |> PEtabODEProblem
stage_problems = PEtabCurriculumProblem(petab_prob, SplitUniform(4))

mdf = petab_prob.model_info.model.petab_tables[:measurements]
unique_t = PEtabTraining._get_unique_timepoints(mdf)
splits = PEtabTraining._makechunks(unique_t, 4; overlap = 1)
model = petab_prob.model_info.model
speciemap = model.speciemap
specie_ids = replace.(string.(first.(speciemap)), ("(t)") => "")
conditions_df, parameters_df = copy(model.petab_tables[:conditions]),
copy(model.petab_tables[:parameters])

# The first window is special, because each parameter already has an initial value
for (i, specie_id) in pairs(specie_ids)
    specie_id in names(conditions_df) && continue
    for (j, cid) in pairs(conditions_df.conditionId)
        # If initial value is an equation expression. In this case NaN must be set to
        # ensure the SBML formula is used for initial value computation
        u0_value = string(speciemap[i].second)
        if !(PEtab.is_number(u0_value) || u0_value in parameters_df.parameterId)
            u0_value = "NaN"
        end
        if specie_id in names(conditions_df)
            conditions_df[j, specie_id] = u0_value
        else
            conditions_df[!, specie_id] .= u0_value
        end
    end
end

observable_df = copy(model.petab_tables[:observables])
measurements_df = copy(model.petab_tables[:measurements])
# Initial values for each window must be added as parameters, window penalty must be
# added as observables + measurement points, and this must be done for each condition
cids = deepcopy(conditions_df.conditionId)
for _cid in cids
    for i in 2:length(splits)
        cid = "__window$(i)_$(_cid)__"
        df_c = conditions_df[conditions_df.conditionId .== _cid, :] |> copy
        df_c[1, :conditionId] = cid
        for specie_id in specie_ids
            pid = "__window$(i)_$(_cid)__$(specie_id)"
            obs_id = "__window$(i)_obs_$(_cid)__$(specie_id)"

            # Set parameters
            df_ps = DataFrame(parameterId = pid, parameterScale = "lin", lowerBound = 0.0,
                upperBound = Inf, nominalValue = 0.0, estimate = 1)
            if "parameterName" in names(parameters_df)
                df_ps[!, :parameterName] .= missing
            end
            parameters_df = vcat(parameters_df, df_ps)

            # Assign initial value to parameter value via conditions table
            df_c[1, specie_id] = pid

            # Setup observables. Multiple window penalty is given by lambda_sqrt as the
            # parameter will be squared during the likelihood computations.
            df_obs = DataFrame(observableId = obs_id,
                observableFormula = "lambda_sqrt * ($(specie_id) - $(obs_id))",
                noiseFormula = "1.0",
                observableTransformation = "lin", noiseDistribution = "normal")
            if "observableName" in names(observable_df)
                df_obs[!, :observableName] .= missing
            end
            observable_df = vcat(observable_df, df_obs)

            # Add in the measurement table the penalty. Given the observable formula above,
            # the measurement value is set to zero, only fulfilled upon merging.
            df_m = DataFrame(observableId = obs_id, simulationConditionId = cid,
                measurement = 0.0, time = maximum(splits[i]))
            if "preequilibrationConditionId" in names(measurements_df)
                df_m[!, :preequilibrationConditionId] .= missing
            end
            if "observableParameters" in names(measurements_df)
                df_m[!, :observableParameters] .= missing
            end
            if "noiseParameters" in names(measurements_df)
                df_m[!, :noiseParameters] .= missing
            end
            if "datasetId" in names(measurements_df)
                df_m[!, :datasetId] .= missing
            end
            measurements_df = vcat(measurements_df, df_m)
        end
        conditions_df = vcat(conditions_df, df_c)
    end
end
df_ps = DataFrame(parameterId = "lambda_sqrt", parameterScale = "lin", lowerBound = 0.0,
    upperBound = Inf, nominalValue = 1.0, estimate = 0)
if "parameterName" in names(parameters_df)
    df_ps[!, :parameterName] .= missing
end
parameters_df = vcat(parameters_df, df_ps)

# Next step. In PEtab.jl add:
#   1. Ability to set t0 for a condition. Will be added to SimulationInfo, and will be manipulated by PEtabTraining.
#   2. Handling of NaN for a non pre-eq. Discussed above, will need to implement an if statement
#   3. Lower level interface of _PEtabModel taking the tables as input. Will remove a lot of code from here.
