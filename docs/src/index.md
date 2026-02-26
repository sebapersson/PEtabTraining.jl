# PEtabTraining.jl

PEtabTraining.jl implements training strategies to make parameter estimation more
efficient for both ordinary differential equation (ODE) models and scientific machine
learning (SciML) models (ODEs with neural network components) more efficient. The package
is designed to be used with [PEtab.jl](https://github.com/sebapersson/PEtab.jl) for ODE
parameter estimation: training strategies can be applied directly to a `PEtabODEProblem` to
obtain modified training objectives (e.g. a multiple-shooting objective) with a single line
of code.

Currently, three training strategies are implemented:

- **Curriculum learning**: strategy where problem difficulty is progressively increased
  across curriculum stages. For ODE problems, this is typically done by gradually
  increasing the number of measurement time points (and often the simulation end time)
  over a fixed number of stages.
- **Multiple shooting**: strategy where the ODE simulation time span is split into windows
  that are fitted jointly. Each window has its own estimated initial state, and a
  continuity penalty is used to promote continuity between adjacent windows.
- **Curriculum multiple shooting**: strategy combining multiple shooting with a curriculum
  schedule. Training starts from a multiple-shooting formulation (often easier to
  optimize) and progressively reduces the number of windows until the original
  single-window problem is recovered.

## Getting started

PEtabTraining.jl is meant to be used together with
[PEtab.jl](https://github.com/sebapersson/PEtab.jl), and all training strategies operate
on a `PEtabODEProblem`. Concrete examples of how to apply the training strategies are
available in the [PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/).
A detailed overview of available options for each training method is provided in the
[API](@ref API) documentation.

## Installation

PEtabTraining.jl is an officially registered Julia package, tested and supported on Linux,
macOS and Windows. The easiest way to install it is via the Julia package manager. In the
Julia REPL, enter:

```julia
julia> ] add PEtabTraining
```

or alternatively

```julia
julia> using Pkg; Pkg.add("PEtabTraining")
```

PEtabTraining is compatible with Julia **1.10 and above**. For best performance, we
strongly recommend using the latest Julia version, which can be most reliably installed
using [juliaup](https://github.com/JuliaLang/juliaup).

## Getting help

If you have any problems using PEtabTraining, here are some helpful tips:

- Post your questions in the `#sciml-sysbio` channel on the
  [Julia Slack](https://julialang.org/slack/).
- If you encounter unexpected behavior or a bug, please see how to file an issue on the
  [Contributing page](https://sebapersson.github.io/PEtabTraining.jl/CONTRIBUTING.md).
