function _PEtabODEProblem(
        model::PEtabModel, prob_original::PEtabODEProblem)::PEtabODEProblem
    @unpack (solver, solver_gradient, ss_solver, ss_solver_gradient, gradient_method,
    hessian_method, sensealg, reuse_sensitivities) = prob_original.probinfo
    prob = PEtabODEProblem(model; odesolver = solver, odesolver_gradient = solver_gradient,
        ss_solver = ss_solver, ss_solver_gradient = ss_solver_gradient,
        gradient_method = gradient_method, hessian_method = hessian_method,
        sensealg = sensealg, reuse_sensitivities = reuse_sensitivities)
    # Sometimes parameters only appear for a subset of splits, for example
    # observable-parameters associated with observables for a subset of conditions and/or
    # late time-points. This can change the order of prob.xnames compared to
    # prob_original.xnames. To avoid problems arising from this, the input for the
    # sub-problems is assumed to follow the order of prob_original, and internally
    # any input x is mapped to the correct order
    prob_to_original = [findfirst(x -> x == name, prob.xnames)
                        for name in prob_original.xnames]
    original_to_prob = [findfirst(x -> x == name, prob_original.xnames)
                        for name in prob.xnames]
    nllh = _get_nllh(prob, original_to_prob)
    prior = _get_prior(prob, original_to_prob)
    chi2 = _get_chi2(prob, original_to_prob)
    simulated_values = _get_simulated_values(prob, original_to_prob)
    residuals = _get_residuals(prob, original_to_prob)
    nllh_grad = _get_nllh_grad(prob, original_to_prob, prob_to_original)
    grad = _get_grad(prob, original_to_prob, prob_to_original)
    grad! = _get_grad!(prob, original_to_prob, prob_to_original)
    grad_prior = _get_grad_prior(prob, original_to_prob, prob_to_original)
    hess = _get_hess(prob, original_to_prob, prob_to_original)
    hess! = _get_hess!(prob, original_to_prob, prob_to_original)
    hess_prior = _get_hess_prior(prob, original_to_prob, prob_to_original)
    return PEtabODEProblem(
        nllh, chi2, grad!, grad, hess!, hess, hess!, hess, nllh_grad, prior,
        grad_prior, hess_prior, simulated_values, residuals, prob.probinfo,
        prob.model_info, prob.nparameters_estimate, prob_original.xnames,
        prob_original.xnominal, prob_original.xnominal_transformed,
        prob_original.lower_bounds, prob_original.upper_bounds)
end

function _get_nllh(prob::PEtabODEProblem, original_to_prob::Vector{Int64})::Function
    nllh = let _original_to_prob = original_to_prob, _prob = prob
        (x) -> begin
            _x = x[_original_to_prob]
            return _prob.nllh(_x)
        end
    end
    return nllh
end

function _get_prior(prob::PEtabODEProblem, original_to_prob::Vector{Int64})::Function
    prior = let _original_to_prob = original_to_prob, _prob = prob
        (x) -> begin
            _x = x[_original_to_prob]
            return _prob.prior(collect(_x))
        end
    end
    return prior
end

function _get_chi2(prob::PEtabODEProblem, original_to_prob::Vector{Int64})::Function
    chi2 = let _original_to_prob = original_to_prob, _prob = prob
        (x) -> begin
            _x = x[_original_to_prob]
            return _prob.chi2(_x)
        end
    end
    return chi2
end

function _get_residuals(prob::PEtabODEProblem, original_to_prob::Vector{Int64})::Function
    residuals = let _original_to_prob = original_to_prob, _prob = prob
        (x) -> begin
            _x = x[_original_to_prob]
            return _prob.residuals(_x)
        end
    end
    return residuals
end

function _get_simulated_values(
        prob::PEtabODEProblem, original_to_prob::Vector{Int64})::Function
    simulated_values = let _original_to_prob = original_to_prob, _prob = prob
        (x) -> begin
            _x = x[_original_to_prob]
            return _prob.simulated_values(_x)
        end
    end
    return simulated_values
end

function _get_grad(prob::PEtabODEProblem, original_to_prob::Vector{Int64},
        prob_to_original::Vector{Int64})::Function
    grad = let _original_to_prob = original_to_prob, _prob_to_original = prob_to_original,
        _prob = prob

        (x) -> begin
            _x = x[_original_to_prob]
            g = _prob.grad(_x)
            return g[_prob_to_original]
        end
    end
    return grad
end

function _get_grad!(prob::PEtabODEProblem, original_to_prob::Vector{Int64},
        prob_to_original::Vector{Int64})::Function
    grad! = let _original_to_prob = original_to_prob, _prob_to_original = prob_to_original,
        _prob = prob

        (g, x) -> begin
            _x = x[_original_to_prob]
            _prob.grad!(g, _x)
            g .= g[_prob_to_original]
            return nothing
        end
    end
    return grad!
end

function _get_nllh_grad(prob::PEtabODEProblem, original_to_prob::Vector{Int64},
        prob_to_original::Vector{Int64})::Function
    nllh_grad = let _original_to_prob = original_to_prob,
        _prob_to_original = prob_to_original, _prob = prob

        (x) -> begin
            _x = x[_original_to_prob]
            nllh, grad = _prob.nllh_grad(_x)
            return nllh, grad[_prob_to_original]
        end
    end
    return nllh_grad
end

function _get_grad_prior(prob::PEtabODEProblem, original_to_prob::Vector{Int64},
        prob_to_original::Vector{Int64})::Function
    grad_prior = let _original_to_prob = original_to_prob,
        _prob_to_original = prob_to_original, _prob = prob

        (g, x) -> begin
            _x = x[_original_to_prob]
            g = _prob.grad_prior(collect(_x))
            g .= g[_prob_to_original]
            return nothing
        end
    end
    return _get_grad_prior
end

function _get_hess(prob::PEtabODEProblem, original_to_prob::Vector{Int64},
        prob_to_original::Vector{Int64})::Function
    hess = let _original_to_prob = original_to_prob, _prob_to_original = prob_to_original,
        _prob = prob

        (x) -> begin
            _x = x[_original_to_prob]
            H = _prob.hess(_x)
            H_out = similar(H)
            _map_hessian!(H_out, H, _prob_to_original)
            return H_out
        end
    end
    return hess
end

function _get_hess!(prob::PEtabODEProblem, original_to_prob::Vector{Int64},
        prob_to_original::Vector{Int64})::Function
    hess! = let _original_to_prob = original_to_prob, _prob_to_original = prob_to_original,
        _prob = prob

        (H_out, x) -> begin
            _x = x[_original_to_prob]
            H = _prob.hess(_x)
            _map_hessian!(H_out, H, _prob_to_original)
            return nothing
        end
    end
    return hess!
end

function _get_hess_prior(prob::PEtabODEProblem, original_to_prob::Vector{Int64},
        prob_to_original::Vector{Int64})::Function
    hess_prior = let _original_to_prob = original_to_prob,
        _prob_to_original = prob_to_original, _prob = prob

        (H_out, x) -> begin
            _x = x[_original_to_prob]
            H_prior = _prob.hess_prior(collect(_x))
            H_out = similar(H_prior)
            _map_hessian!(H_out, H_prior, _prob_to_original)
            return H_out
        end
    end
    return hess_prior
end

function _map_hessian!(H_original::T, H_prob::T,
        prob_to_original::Vector{Int64})::Nothing where {T <: Matrix{<:AbstractFloat}}
    for (i1, i2) in pairs(prob_to_original)
        for (j1, j2) in pairs(prob_to_original)
            H_original[i1, j1] = H_prob[i2, j2]
        end
    end
    return nothing
end
