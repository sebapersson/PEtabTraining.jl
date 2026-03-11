function PEtabClMsProblem(
        prob_original::PEtabODEProblem, split_alg::SplitTime;
        regularization_obs::Union{Nothing, String, Symbol} = nothing,
        regularization_specie::Union{Nothing, String, Symbol} = nothing,
        window_u0_scale::Symbol = :lin, window_penalty_scale::Symbol = :lin
    )::PEtabClMsProblem
    if prob_original.model_info.simulation_info.has_pre_equilibration
        throw(ArgumentError("Curriculum + multiple shooting approach is not supported for \
            models with pre-equilibration simulation conditions."))
    end

    @argcheck window_u0_scale in [:lin, :log, :log10]
    @argcheck window_penalty_scale in [:lin, :log, :log10]

    _check_regularization_specie(regularization_obs, regularization_specie)

    ms_windows_stages = Dict{Symbol, Vector{Vector{Float64}}}()
    ms_windows_stages[:stage1] = _split_ms(split_alg, prob_original)

    n_windows = length(ms_windows_stages[:stage1])
    for i in 2:n_windows
        ms_windows_prev = ms_windows_stages[Symbol("stage$(i - 1)")]
        ms_windows_stage = Vector{Vector{Float64}}(undef, n_windows - i + 1)
        for j in eachindex(ms_windows_stage)
            ms_windows_stage[j] = unique(reduce(vcat, ms_windows_prev[j:(j + 1)]))
            ms_windows_stage[j] = [first(ms_windows_stage[j]), last(ms_windows_stage[j])]
        end
        ms_windows_stages[Symbol("stage$(i)")] = ms_windows_stage
    end

    petab_problems = Vector{PEtabODEProblem}(undef, n_windows)
    for i in 1:(n_windows - 1)
        ms_windows_stage = ms_windows_stages[Symbol("stage$(i)")]
        petab_problems[i] = PEtabTraining._get_petab_ms_problem(
            prob_original, ms_windows_stage, window_u0_scale, window_penalty_scale,
            _string(regularization_obs), _string(regularization_specie)
        )
    end
    petab_problems[n_windows] = prob_original
    return PEtabClMsProblem(
        petab_problems, prob_original, split_alg, ms_windows_stages, window_u0_scale,
        window_penalty_scale
    )
end
