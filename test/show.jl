using PEtabTraining, Test

include(joinpath(@__DIR__, "common.jl"))

# Basic show for each problem type
model_id = "Boehm_JProteomeRes2014"
prob_original = _get_petab_problem(model_id)
prob_cl = PEtabClProblem(prob_original, SplitTime(3))
@test "$prob_cl" == "PEtabClProblem Boehm_JProteomeRes2014: 3 CL stages, 9 parameters\n(for \
    more statistics, call `describe(prob_cl)`)"
prob_ms = PEtabMsProblem(prob_original, SplitTime(3))
@test "$prob_ms" == "PEtabMsProblem Boehm_JProteomeRes2014: 3 MS windows, 9 \
    parameters\n(for more statistics, call `describe(prob_ms)`)"
prob_cl_ms = PEtabClMsProblem(prob_original, SplitTime(3))
@test "$prob_cl_ms" == "PEtabClMsProblem Boehm_JProteomeRes2014: 3 CL stages, 9 \
    parameters\n(for more statistics, call `describe(prob_cl_ms)`)"

# Describe for models with and without simulation conditions
# Curriculum
model_id = "Boehm_JProteomeRes2014"
prob_original = _get_petab_problem(model_id)
prob_cl1 = PEtabClProblem(prob_original, SplitTime(3))
model_id = "Bachmann_MSB2011"
prob_original = _get_petab_problem(model_id)
prob_cl2 = PEtabClProblem(prob_original, SplitTime(4))
@test "$(describe(prob_cl1; as_string = true))" == "PEtabClProblem Boehm_JProteomeRes2014\
    \nProblem statistics\n  Parameters to estimate: 9\n  ODE: 8 states, 10 parameters\n  \
    Observables: 3\n  Simulation conditions: 1\n\nCurriculum statistics (3 stages)\n  \
    Stage 1: tspan [0.0, 20.0]\n           fraction (obs/cond): 1.0/1.0\n  \
    Stage 2: tspan [0.0, 80.0]\n           fraction (obs/cond): 1.0/1.0\n  \
    Stage 3: tspan [0.0, 240.0]\n           fraction (obs/cond): 1.0/1.0"
@test "$(describe(prob_cl2; as_string = true))" == "PEtabClProblem Bachmann_MSB2011\
    \nProblem statistics\n  Parameters to estimate: 113\n  ODE: 25 states, 39 parameters\n  \
    Observables: 20\n  Simulation conditions: 36\n\nCurriculum statistics (4 stages)\n  \
    Stage 1: tspan [0.0, 7.5]\n           fraction (obs/cond): 0.5/0.5\n  \
    Stage 2: tspan [0.0, 25.0]\n           fraction (obs/cond): 0.69/0.8\n  \
    Stage 3: tspan [0.0, 80.0]\n           fraction (obs/cond): 0.86/0.8\n  \
    Stage 4: tspan [0.0, 360.0]\n           fraction (obs/cond): 1.0/1.0"
@test "$(describe(prob_ms; as_string = true))" == "PEtabMsProblem Boehm_JProteomeRes2014\
    \nProblem statistics\n  Parameters to estimate: 9\n  ODE: 8 states, 10 parameters\n  \
    Observables: 3\n  Simulation conditions: 1\n\nWindow statistics (3 windows)\n  \
    Penalty λ = 1.0e+00, window u0 parameters: 16\n  Window 1: tspan [0.0, 30.0]\n  \
    Window 2: tspan [30.0, 100.0]\n  Window 3: tspan [100.0, 240.0]"
@test "$(describe(prob_cl_ms; as_string = true))" == "PEtabClMsProblem \
    Boehm_JProteomeRes2014\nProblem statistics\n  Parameters to estimate: 9\n  ODE: 8 \
    states, 10 parameters\n  Observables: 3\n  Simulation conditions: 1\n\nCurriculum \
    statistics (3 stages)\n  Window penalty λ = 1.0e+00\n  Stage 1: window tspans \
    [0, 30], [30, 100], [100, 240]\n  Stage 2: window tspans [0, 100], [30, 240]\n  \
    Stage 3: original problem"
