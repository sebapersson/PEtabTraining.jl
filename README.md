# PEtabTraining.jl
*Training strategies for parameter estimation in dynamic (ODE) models*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sebapersson.github.io/PEtabTraining.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sebapersson.github.io/PEtabTraining.jl/dev/)
[![Build Status](https://github.com/sebapersson/PEtabTraining.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sebapersson/PEtabTraining.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/sebapersson/PEtabTraining.jl/graph/badge.svg?token=1J7RKUPZOM)](https://codecov.io/gh/sebapersson/PEtabTraining.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

[Documentation](https://sebapersson.github.io/PEtabTraining.jl/stable/) |
[Contributing](https://sebapersson.github.io/PEtabTraining.jl/CONTRIBUTING.md)

PEtabTraining.jl implements training strategies to make parameter estimation more efficient
for both ordinary differential equation (ODE) models and scientific machine learning
(SciML) models (ODEs with neural network components) more efficient. The package is
designed to be used with [PEtab.jl](https://github.com/sebapersson/PEtab.jl) for ODE
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

Concrete examples on how to apply these training strategies can be found in the
[PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/), while the
[PEtabTraining.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/) holds
the detailed options for each strategy.

## Installation

PEtabTraining.jl is a registered Julia package and can be installed via the Julia package
manager:

```julia
import Pkg
Pkg.add("PEtabTraining")
```

PEtabTraining.jl is compatible with Julia 1.10 and above. For additional installation
details, see the online
[documentation](https://sebapersson.github.io/PEtabTraining.jl/stable/).
