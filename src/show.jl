import Base

function Base.show(io::Core.IO, prob_cl::PEtabClProblem)
    nx = length(prob_cl.original.xnominal)
    name = prob_cl.original.model_info.model.name
    n_stages = length(prob_cl.petab_problems)

    header = styled"{bold:PEtabClProblem} {emphasis:$(name)}: $(n_stages) CL stages, $nx \
        parameters\n(for more statistics, call `describe(prob_cl)`)"
    return print(io, styled"$(header)")
end
function Base.show(io::Core.IO, prob_ms::PEtabMsProblem)
    nx = length(prob_ms.original.xnominal)
    name = prob_ms.original.model_info.model.name
    n_windows = length(prob_ms.ms_windows)

    header = styled"{bold:PEtabMsProblem} {emphasis:$(name)}: $(n_windows) MS windows, $nx \
        parameters\n(for more statistics, call `describe(prob_ms)`)"
    return print(io, styled"$(header)")
end
function Base.show(io::Core.IO, prob_cl_ms::PEtabClMsProblem)
    nx = length(prob_cl_ms.original.xnominal)
    name = prob_cl_ms.original.model_info.model.name
    n_stages = length(prob_cl_ms.petab_problems)

    header = styled"{bold:PEtabClMsProblem} {emphasis:$(name)}: $(n_stages) CL stages, $nx \
        parameters\n(for more statistics, call `describe(prob_cl_ms)`)"
    return print(io, styled"$(header)")
end

"""
    describe(prob::PEtabClProblem)

Print summary statistics for the overall estimation problem and for each curriculum stage.

Fractions are reported relative to the original problem (fraction of observables / simulation
conditions included in each stage).
"""
function StatsBase.describe(prob_cl::PEtabClProblem; as_string::Bool = false)
    name = prob_cl.original.model_info.model.name
    header = styled"{bold:PEtabClProblem} {emphasis:$(name)}\n"
    model_stats = _get_model_stats(prob_cl.original)
    model_describe = _describe_common(prob_cl.original)

    # Stats for each CL-stage
    n_stages = length(prob_cl.petab_problems)
    cl_stats = styled"{underline:Curriculum statistics} ($(n_stages) stages)\n"
    for (i, prob_stage) in pairs(prob_cl.petab_problems)
        stage_stats = _get_model_stats(prob_stage)
        time_span = round.(_get_time_span(prob_stage), digits = 2)

        frac_conditions = round(
            stage_stats.n_conditions / model_stats.n_conditions; digits = 2
        )
        frac_observables = round(
            stage_stats.n_observables / model_stats.n_observables; digits = 2
        )

        cl_stats *= "  Stage $i: tspan [$(time_span[1]), $(time_span[2])]\n"
        cl_stats *= "           fraction (obs/cond): $(frac_conditions)/$(frac_observables)"
        if i != length(prob_cl.petab_problems)
            cl_stats *= "\n"
        end
    end

    if as_string
        return "$(header)$(model_describe)$(cl_stats)"
    end
    return print(styled"$(header)$(model_describe)$(cl_stats)")
end
"""
    describe(prob::PEtabMsProblem)

Print summary statistics for the multiple-shooting estimation problem and window
configuration (penalty, time spans, and number of initial window parameters).
"""
function StatsBase.describe(prob_ms::PEtabMsProblem; as_string::Bool = false)
    @unpack original, ms_windows, petab_ms_problem = prob_ms

    name = original.model_info.model.name
    header = styled"{bold:PEtabMsProblem} {emphasis:$(name)}\n"
    model_describe = _describe_common(original)

    n_windows = length(ms_windows)
    n_initial_parameters = length(_get_ms_u0_x_names(petab_ms_problem))
    window_penalty = @sprintf("%.1e", _get_ms_window_penalty(petab_ms_problem))

    ms_stats = styled"{underline:Window statistics} ($(n_windows) windows)\n"
    ms_stats *= "  Penalty λ = $(window_penalty), window u0 parameters: \
        $(n_initial_parameters)\n"
    for (i, ms_window) in pairs(ms_windows)
        ms_stats *= "  Window $i: tspan [$(ms_window[1]), $(ms_window[2])]"
        if i != n_windows
            ms_stats *= '\n'
        end
    end

    if as_string
        return "$(header)$(model_describe)$(ms_stats)"
    end
    return print(styled"$(header)$(model_describe)$(ms_stats)")
end
"""
    describe(prob::PEtabClMsProblem)

Print summary statistics for the overall estimation problem and for each curriculum stage.

For each stage, multiple-shooting window statistics are shown (window time-span and window
penalty λ).
"""
function StatsBase.describe(prob_cl_ms::PEtabClMsProblem; as_string::Bool = false)
    name = prob_cl_ms.original.model_info.model.name
    header = styled"{bold:PEtabClMsProblem} {emphasis:$(name)}\n"
    model_stats = _get_model_stats(prob_cl_ms.original)
    model_describe = _describe_common(prob_cl_ms.original)

    # Stats for each CL-stage
    window_penalty = @sprintf("%.1e", _get_ms_window_penalty(prob_cl_ms.petab_problems[1]))
    n_stages = length(prob_cl_ms.petab_problems)
    cl_stats = styled"{underline:Curriculum statistics} ($(n_stages) stages)\n"
    cl_stats *= "  Window penalty λ = $(window_penalty)\n"
    for i in eachindex(prob_cl_ms.petab_problems)
        if i == length(prob_cl_ms.petab_problems)
            cl_stats *= "  Stage $i: original problem"
            continue
        end

        ms_windows = prob_cl_ms.ms_windows[Symbol("stage$i")]
        ms_windows_fmt = _fmt_windows(ms_windows)
        cl_stats *= "  Stage $i: window tspans $(ms_windows_fmt)\n"
    end
    if as_string
        return "$(header)$(model_describe)$(cl_stats)"
    end
    return print(styled"$(header)$(model_describe)$(cl_stats)")
end

function _describe_common(prob::PEtabODEProblem)
    model_stats = _get_model_stats(prob)

    opt_head = styled"{underline:Problem statistics}\n"
    opt1 = "  Parameters to estimate: $(model_stats.nx)\n"
    opt2 = "  ODE: $(model_stats.n_states) states, $(model_stats.n_ps_ode) parameters\n"
    opt3 = "  Observables: $(model_stats.n_observables)\n"
    opt4 = "  Simulation conditions: $(model_stats.n_conditions)\n"
    return styled"$(opt_head)$(opt1)$(opt2)$(opt3)$(opt4)\n"
end

function _get_model_stats(prob::PEtabODEProblem)
    model = prob.model_info.model
    nx = length(prob.xnominal)
    n_states = length(PEtab._get_state_ids(model.sys_mutated))
    n_ps_ode = PEtab._get_n_parameters_sys(model.sys_mutated)
    n_observables = length(unique(model.petab_tables[:measurements].observableId))
    n_conditions = length(prob.model_info.simulation_info.conditionids[:experiment])

    return (
        nx = nx, n_states = n_states, n_ps_ode = n_ps_ode, n_observables = n_observables,
        n_conditions = n_conditions,
    )
end

function _get_time_span(prob::PEtabODEProblem)
    measurements_df = prob.model_info.model.petab_tables[:measurements]
    return [minimum(measurements_df.time), maximum(measurements_df.time)]
end

_fmt_num(x) = @sprintf("%.3g", float(x))
_fmt_span(t0, t1) = "[" * _fmt_num(t0) * ", " * _fmt_num(t1) * "]"

function _fmt_windows(ms_windows::Vector{<:Vector{<:Real}}; max_show::Int = 5)
    n = length(ms_windows)
    spans = map(w -> _fmt_span(w[1], w[2]), ms_windows)

    if n ≤ max_show
        return join(spans, ", ")
    else
        # show first 2 and last 2
        head = spans[1:2]
        tail = spans[(end - 1):end]
        return join(vcat(head, ["…"], tail), ", ")
    end
end
