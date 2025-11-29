"""
    set_u0_windows!(prob::Union{PEtabMSProblem, PEtabMSProblem}, method::Symbol, value)
"""
function set_u0_windows!(prob::PEtabCLMSProblem, method::Symbol, x)::Nothing
    @unpack petab_problems, windows, original = prob
    for i in 1:(length(petab_problems) - 1)
        windows_stage = windows[Symbol("stage$(i)")]
        _set_u0_windows!(petab_problems[i], original, windows_stage, method, x)
    end
    return nothing
end
function set_u0_windows!(prob::PEtabMSProblem, method::Symbol, x)::Nothing
    @unpack petab_prob_ms, original, windows = prob
    _set_u0_windows!(petab_prob_ms, original, windows, method, x)
    return nothing
end

function _set_u0_windows!(prob::PEtabODEProblem, prob_original::PEtabODEProblem,
        windows::Vector{Vector{Float64}}, method::Symbol, x)::Nothing
    @argcheck method in [:constant, :window1_u0, :window1_simulate]
    if method == :constant
        @argcheck isa(x, Real) "For method :constant, x must be a Real value"
    else
        @argcheck isa(x, ComponentArrays.ComponentVector) || isa(x, AbstractVector) "For \
            method :window1_*, x must be a ComponentVector or a Vector"
    end
    if method == :constant
        _set_u0_windows_as_constant!(prob, x)
    elseif method == :window1_u0
        _set_u0_windows_as_window1(prob, prob_original, x)
    else
        _set_u0_windows_as_simulate!(prob, prob_original, windows, x)
    end
    _transform_x!(prob)
    return nothing
end

function _set_u0_windows_as_constant!(prob::PEtabODEProblem, value::Real)::Nothing
    u0_xnames = _get_ms_u0_xnames(prob)
    @views prob.xnominal[u0_xnames] .= value
    return nothing
end

function _set_u0_windows_as_window1(prob::PEtabODEProblem, prob_original::PEtabODEProblem, x)::Nothing
    specie_ids = _get_specie_ids(prob_original)
    xnames_u0 = _get_ms_u0_xnames(prob)
    for xname in xnames_u0
        cid = _get_cid_from_window_id(xname)
        specie_id = _get_specie_id_from_window_id(xname)
        u0_cid = PEtab.get_u0(x, prob_original; retmap = false, cid = cid)
        u0_value = u0_cid[findfirst(x -> x == specie_id, specie_ids)]
        @views prob.xnominal[xname] = u0_value
    end
    return nothing
end

function _set_u0_windows_as_simulate!(
        prob::PEtabODEProblem, prob_original::PEtabODEProblem,
        windows::Vector{Vector{Float64}}, x)::Nothing
    specie_ids = _get_specie_ids(prob_original)
    xnames_u0 = _get_ms_u0_xnames(prob)
    cids_original = string.(prob_original.model_info.simulation_info.conditionids[:experiment])
    for cid in cids_original
        sol = PEtab.get_odesol(x, prob_original; cid = cid)
        for xname in xnames_u0
            _cid = _get_cid_from_window_id(xname)
            _cid != cid && continue
            specie_id = _get_specie_id_from_window_id(xname)
            specie_index = findfirst(x -> x == specie_id, specie_ids)
            t0_window = minimum(windows[_get_index_from_window_id(xname)])
            @views prob.xnominal[xname] = sol(t0_window)[specie_index]
        end
    end
    return nothing
end

function set_window_penalty!(prob::PEtabCLMSProblem, x::Real)::Nothing
    for i in 1:(length(prob.petab_problems) - 1)
        _set_window_penalty!(prob.petab_problems[i], x)
    end
    return nothing
end
function set_window_penalty!(prob::PEtabMSProblem, x::Real)::Nothing
    _set_window_penalty!(prob.petab_prob_ms, x)
    return nothing
end

function _set_window_penalty!(prob::PEtabODEProblem, x::Real)::Nothing
    @argcheck x≥0 "Multiple shooting window penalty parameter must be ≥0"
    petab_parameters = prob.model_info.petab_parameters
    ix = findfirst(x -> x == :lambda_sqrt, petab_parameters.parameter_id)
    petab_parameters.nominal_value[ix] = sqrt(x)
    return nothing
end

"""
map_x_stage(x, clms::PEtabCLMSProblem; from::Integer=1, to::Integer=2)

Map input Vector `x` coming from stage `from` to the layout of stage `to` for a
`PEtabCLMSProblem`

`from` and `to` must be valid stage indices. `x` can be a `Vector` or a `ComponentVector`,
with the ordering expected by `PEtabODEproblem` in stage `from`.
"""
function map_x_stage(x, clms::PEtabCLMSProblem, from::Integer = 1, to::Integer = 2)
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
