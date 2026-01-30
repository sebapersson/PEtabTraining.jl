"""
    set_u0_ms_windows!(x, prob_ms::PEtabMsProblem; init = MsInitConstant(0.01))
    set_u0_ms_windows!(x, prob_ms::PEtabMsProblem, p; init = MsInitFirst())
    set_u0_ms_windows!(x, prob_ms::PEtabMsProblem, p; init = MsInitSimulate())

Set multiple-shooting window initial values in the multiple-shooting parameter vector `x`.

`x` must be a parameter vector in the format expected by `prob_ms` (e.g. returned by
`get_x(prob_ms)`). When provided, `p` must be a parameter vector compatible with either the
original PEtab problem or the multiple-shooting PEtab problem.

Initialization is controlled by `init`:
- `MsInitConstant(value)`: Set all window initial values to `value`.
- `MsInitFirst()`: Copy initial values from the first window to all windows.
- `MsInitSimulate(...)`: Forward simulate the model and set each window initial value
  from the simulated state at the window start time.
"""
function set_u0_ms_windows!(
        x::Union{Vector{<:Real}, ComponentArray}, prob::PEtabMsProblem;
        init::MsInitConstant = MsInitConstant()
    )::Nothing
    @unpack petab_ms_problem = prob
    _set_u0_ms_windows!(x, petab_ms_problem, init)
    return nothing
end
function set_u0_ms_windows!(
        x::Union{Vector{<:Real}, ComponentArray}, prob::PEtabMsProblem,
        p::Union{Vector{<:Real}, ComponentArray};
        init::Union{MsInitFirst, MsInitSimulate} = MsInitFirst()
    )::Nothing

    @unpack petab_ms_problem, original, ms_windows = prob
    return _set_u0_ms_windows!(x, petab_ms_problem, original, p, ms_windows, init)
end
"""
    set_u0_ms_windows!(x, prob::PEtabClMSProblem, stage::Integer; init = MsInitConstant(0.01))
    set_u0_ms_windows!(x, prob::PEtabClMSProblem, stage::Integer, p; init = MsInitFirst())
    set_u0_ms_windows!(x, prob::PEtabClMSProblem, stage::Integer, p; init = MsInitSimulate())

Set multiple-shooting window initial values in the multiple-shooting parameter vector `x`
for `stage` of a combined curriculum-learning and multiple-shooting problem `prob`.

`x` must be in the format expected by the multiple-shooting stage problem (e.g. returned by
`get_x(prob; stage = stage)`). All other arguments have the same meaning as for
`set_u0_ms_windows!(x, prob::PEtabMsProblem, ...)`.
"""
function set_u0_ms_windows!(
        x::Union{Vector{<:Real}, ComponentArray}, prob::PEtabClMsProblem, stage::Integer;
        init::MsInitConstant = MsInitConstant()
    )::Nothing
    @argcheck 1 ≤ stage ≤ length(prob.petab_problems) "invalid `stage` index"

    if stage == length(prob.petab_problems)
        return nothing
    end

    petab_ms_problem = prob.petab_problems[stage]
    _set_u0_ms_windows!(x, petab_ms_problem, init)
    return nothing
end
function set_u0_ms_windows!(
        x::Union{Vector{<:Real}, ComponentArray}, prob::PEtabClMsProblem, stage::Integer,
        p::Union{Vector{<:Real}, ComponentArray};
        init::Union{MsInitFirst, MsInitSimulate} = MsInitFirst()
    )::Nothing
    @argcheck 1 ≤ stage ≤ length(prob.petab_problems) "invalid `stage` index"

    petab_ms_problem = prob.petab_problems[stage]
    ms_windows = prob.ms_windows[Symbol("stage$(stage)")]
    return _set_u0_ms_windows!(x, petab_ms_problem, prob.original, p, ms_windows, init)
end

function _set_u0_ms_windows!(
        x, petab_ms_problem::PEtabODEProblem, init::MsInitConstant
    )::Nothing
    x_names = ComponentArrays.labels(petab_ms_problem.xnominal)

    u0_x_names = _get_ms_u0_x_names(petab_ms_problem)
    idx = [findfirst(x -> x == name, x_names) for name in string.(u0_x_names)]
    u0_values = fill(init.value, length(idx))
    PEtab.transform_x!(u0_values, u0_x_names, petab_ms_problem.model_info.xindices)
    x[idx] .= u0_values
    return nothing
end
function _set_u0_ms_windows!(
        x, petab_ms_problem::PEtabODEProblem, original::PEtabODEProblem, p,
        ::Vector{<:Vector{<:Real}}, ::MsInitFirst
    )::Nothing
    p_original = _get_p_original(p, petab_ms_problem, original)

    specie_ids = _get_specie_ids(original)

    x_names = ComponentArrays.labels(petab_ms_problem.xnominal)
    u0_x_names = _get_ms_u0_x_names(petab_ms_problem)
    for u0_x_name in string.(u0_x_names)
        condition_id = _get_condition_id_from_window(u0_x_name)
        specie_id = _get_specie_id_from_window(u0_x_name)

        u0_condition_id = PEtab.get_u0(
            p_original, original; retmap = false, condition = condition_id
        )
        u0_value = u0_condition_id[findfirst(x -> x == specie_id, specie_ids)]

        idx = findfirst(x -> x == u0_x_name, x_names)
        x[idx] = u0_value
    end
    return nothing
end
function _set_u0_ms_windows!(
        x, petab_ms_problem::PEtabODEProblem, original::PEtabODEProblem, p,
        ms_windows::Vector{<:Vector{<:Real}}, ::MsInitSimulate
    )::Nothing
    p_original = _get_p_original(p, petab_ms_problem, original)

    specie_ids = _get_specie_ids(original)

    x_names = ComponentArrays.labels(petab_ms_problem.xnominal)
    u0_x_names = _get_ms_u0_x_names(petab_ms_problem)

    condition_ids_original = string.(
        original.model_info.simulation_info.conditionids[:experiment]
    )
    for condition_id in condition_ids_original
        sol = PEtab.get_odesol(p_original, original; condition = condition_id)
        for u0_x_name in string.(u0_x_names)
            _condition_id = _get_condition_id_from_window(u0_x_name)
            _condition_id != condition_id && continue

            specie_id = _get_specie_id_from_window(u0_x_name)
            specie_index = findfirst(x -> x == specie_id, specie_ids)
            t0_window = minimum(ms_windows[_get_index_from_window(u0_x_name)])

            idx = findfirst(x -> x == u0_x_name, x_names)
            @views x[idx] = sol(t0_window)[specie_index]
        end
    end
    return nothing
end

function set_window_penalty!(prob::PEtabClMsProblem, x::Real)::Nothing
    for i in 1:(length(prob.petab_problems) - 1)
        _set_window_penalty!(prob.petab_problems[i], x)
    end
    return nothing
end
function set_window_penalty!(prob::PEtabMsProblem, x::Real)::Nothing
    _set_window_penalty!(prob.petab_ms_problem, x)
    return nothing
end

function _set_window_penalty!(prob::PEtabODEProblem, x::Real)::Nothing
    @argcheck x ≥ 0 "Multiple shooting window penalty parameter must be ≥0"
    petab_parameters = prob.model_info.petab_parameters
    ix = findfirst(x -> x == :lambda_sqrt, petab_parameters.parameter_id)
    petab_parameters.nominal_value[ix] = sqrt(x)
    return nothing
end

function _get_p_original(
        p, petab_ms_problem::PEtabODEProblem, original::PEtabODEProblem
    )
    x_original = original.xnominal_transformed
    x_ms = petab_ms_problem.xnominal_transformed

    if length(p) != length(x_original) && length(p) != x_ms
        throw(ArgumentError("Input p must be a Vector/ComponentArray with equal \
            length to either the parameter vector in the original PEtabODEProblem \
            or the PEtabMs/PEtabClMS-Problem"))
    end

    if length(p) == length(x_original)
        return p
    end

    idx = _perm_from_labels(xnominal_transformed, xnominal_transformed)
    return p[idx]
end

"""
map_x_stage(x, clms::PEtabClMsProblem; from::Integer=1, to::Integer=2)

Map input Vector `x` coming from stage `from` to the layout of stage `to` for a
`PEtabClMsProblem`

`from` and `to` must be valid stage indices. `x` can be a `Vector` or a `ComponentVector`,
with the ordering expected by `PEtabODEProblem` in stage `from`.
"""
function map_x_stage(x, clms::PEtabClMsProblem, from::Integer = 1, to::Integer = 2)
    @argcheck 1 ≤ from ≤ length(clms.petab_problems) "invalid `from` stage index"
    @argcheck 1 ≤ to ≤ length(clms.petab_problems) "invalid `to` stage index"
    @argcheck from != to

    # target prototype to preserve types/ordering/defaults
    x_to = clms.petab_problems[to].xnominal_transformed |> deepcopy
    x_from = clms.petab_problems[from].xnominal_transformed |> deepcopy
    if from < to
        ix = _perm_from_labels(x_to, x_from)
        x_to .= x[ix]
    else
        ix = Int64[]
        for label in ComponentArrays.labels(x_from)
            if label in ComponentArrays.labels(x_to)
                push!(ix, only(ComponentArrays.label2index(x_to, label)))
            end
        end
        x_to[ix] .= x
    end
    return x_to
end
