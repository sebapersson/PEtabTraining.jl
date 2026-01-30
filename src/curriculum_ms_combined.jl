struct PEtabCLMSProblem
    split_alg::Any
    windows::Dict{Symbol, Vector{Vector{Float64}}}
    petab_problems::Vector{PEtabODEProblem}
    original::PEtabODEProblem
end
function PEtabCLMSProblem(prob_original::PEtabODEProblem, split_alg;
        regularization_obs::Union{Nothing, String, Symbol} = nothing,
        regularization_specie::Union{Nothing, String, Symbol} = nothing)::PEtabCLMSProblem
    _check_regularization_specie(regularization_obs, regularization_specie)
    windows_stages = Dict{Symbol, Vector{Vector{Float64}}}()
    windows_stages[:stage1] = _split(split_alg, prob_original, :multiple_shooting)

    n_windows = length(windows_stages[:stage1])
    for i in 2:n_windows
        windows_stage_prev = windows_stages[Symbol("stage$(i-1)")]
        windows_stage = Vector{Vector{Float64}}(undef, n_windows - i + 1)
        for j in eachindex(windows_stage)
            windows_stage[j] = unique(reduce(vcat, windows_stage_prev[j:(j + 1)]))
        end
        windows_stages[Symbol("stage$(i)")] = windows_stage
    end

    petab_problems = Vector{PEtabODEProblem}(undef, n_windows)
    for i in 1:(n_windows - 1)
        windows_stage = windows_stages[Symbol("stage$(i)")]
        petab_problems[i] = PEtabTraining._get_petab_ms_problem(
            prob_original, windows_stage, _string(regularization_obs),
            _string(regularization_specie))
    end
    petab_problems[n_windows] = prob_original
    return PEtabCLMSProblem(split_alg, windows_stages, petab_problems, prob_original)
end
