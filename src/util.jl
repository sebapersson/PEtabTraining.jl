"""
    set_u0_windows!(prob::PEtabMultipleShootingProblem, value::Real)::Nothing

Set u0 for each specie for windows the provided value
"""
function set_u0_windows!(prob::PEtabMultipleShootingProblem, value::Real)::Nothing
    xnames_u0 = _get_ms_u0_xnames(prob)
    @views prob.petab_prob_ms.xnominal[xnames_u0] .= value
    _transform_x!(prob)
    return nothing
end
function set_u0_windows!(prob::PEtabMultipleShootingProblem, x, method::Symbol)::Nothing
    @assert method in [:window1_u0, :window1_simulate]
    if method == :window1_u0
        _set_u0_windows_window1_u0(prob, x)
    else
        _set_u0_windows_window1_simulate(prob, x)
    end
    _transform_x!(prob)
    return nothing
end

function set_window_penalty!(prob::PEtabMultipleShootingProblem, x::Real)::Nothing
    @argcheck x≥0 "Multiple shooting window penalty parameter must be ≥0"
    petab_parameters = prob.petab_prob_ms.model_info.petab_parameters
    ix = findfirst(x -> x == :lambda_sqrt, petab_parameters.parameter_id)
    petab_parameters.nominal_value[ix] = sqrt(x)
    return nothing
end

function _set_u0_windows_window1_u0(prob::PEtabMultipleShootingProblem, x)::Nothing
    @unpack original, petab_prob_ms = prob
    specie_ids = _get_specie_ids(original)
    xnames_u0 = _get_ms_u0_xnames(prob)
    for xname in xnames_u0
        cid = _get_cid_from_window_id(xname)
        specie_id = _get_specie_id_from_window_id(xname)
        u0_cid = PEtab.get_u0(x, original; retmap = false, cid = cid)
        u0_value = u0_cid[findfirst(x -> x == specie_id, specie_ids)]
        @views petab_prob_ms.xnominal[xname] = u0_value
    end
    return nothing
end

function _set_u0_windows_window1_simulate(prob::PEtabMultipleShootingProblem, x)::Nothing
    @unpack original, petab_prob_ms = prob
    windows = _split(prob.split_algorithm, original, :multiple_shooting)
    specie_ids = _get_specie_ids(original)
    xnames_u0 = _get_ms_u0_xnames(prob)
    cids_original = string.(original.model_info.simulation_info.conditionids[:experiment])
    for cid in cids_original
        sol = PEtab.get_odesol(x, original; cid = cid)
        for xname in xnames_u0
            _cid = _get_cid_from_window_id(xname)
            _cid != cid && continue
            specie_id = _get_specie_id_from_window_id(xname)
            specie_index = findfirst(x -> x == specie_id, specie_ids)
            t0_window = maximum(windows[_get_index_from_window_id(xname)])
            @views petab_prob_ms.xnominal[xname] = sol(t0_window)[specie_index]
        end
    end
    return nothing
end

function _transform_x!(prob::PEtabMultipleShootingProblem)::Nothing
    @unpack xnominal, xnominal_transformed, xnames, model_info = prob.petab_prob_ms
    @views xnominal_transformed .= PEtab.transform_x(xnominal, xnames, model_info.xindices;
        to_xscale = true)
    return nothing
end

function _get_ms_u0_xnames(prob::PEtabMultipleShootingProblem)::Vector{Symbol}
    xnames = Symbol[]
    for xname in prob.petab_prob_ms.xnames
        xname in prob.original.xnames && continue
        push!(xnames, xname)
    end
    return xnames
end
