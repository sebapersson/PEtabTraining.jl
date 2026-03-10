# PEtabTraining.jl

_Training strategies for parameter estimation in dynamic models_

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sebapersson.github.io/PEtabTraining.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sebapersson.github.io/PEtabTraining.jl/dev/)
[![Build Status](https://github.com/sebapersson/PEtabTraining.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sebapersson/PEtabTraining.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/sebapersson/PEtabTraining.jl/graph/badge.svg?token=1J7RKUPZOM)](https://codecov.io/gh/sebapersson/PEtabTraining.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

[Documentation](https://sebapersson.github.io/PEtabTraining.jl/stable/) |
[Contributing](https://sebapersson.github.io/PEtabTraining.jl/CONTRIBUTING.md)

PEtabTraining.jl implements training strategies that improve the efficiency of parameter
estimation for both ordinary differential equation (ODE) models and scientific machine
learning (SciML) models. The package is designed to be used with
[PEtab.jl](https://github.com/sebapersson/PEtab.jl): training strategies can be applied
directly to a `PEtabODEProblem` to obtain modified training objectives (e.g. multiple
shooting objective) with a single line of code.

Currently, three training strategies are implemented:

- **Curriculum learning**: strategy where problem difficulty is progressively increased
  across curriculum stages. For dynamic models, this is typically done by gradually
  increasing the number of measurement time points (and often the simulation end time) over
  a fixed number of stages.
- **Multiple shooting**: strategy where the ODE simulation time span is split into windows
  that are fitted jointly. Each window has its own estimated initial state, and a continuity
  penalty is used to promote continuity between adjacent windows.
- **Curriculum multiple shooting**: strategy combining multiple shooting with a curriculum
  schedule. Training starts from a multiple-shooting formulation (often easier to optimize)
  and progressively reduces the number of windows until the original single-window problem
  is recovered.

Concrete examples of how to apply these training strategies are available in the
[PEtab.jl documentation](https://sebapersson.github.io/PEtab.jl/stable/), while the
[PEtabTraining.jl documentation](https://sebapersson.github.io/PEtabTraining.jl/stable/)
contains the detailed API and options for each strategy.

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

## Citation

If you use PEtabTraining.jl in work that is published, please cite the paper below:

```bibtex
@article{PEtabBioinformatics2025,
  title={PEtab.jl: advancing the efficiency and utility of dynamic modelling},
  author={Persson, Sebastian and Fr{\"o}hlich, Fabian and Grein, Stephan and Loman, Torkel and Ognissanti, Damiano and Hasselgren, Viktor and Hasenauer, Jan and Cvijovic, Marija},
  journal={Bioinformatics},
  volume={41},
  number={9},
  pages={btaf497},
  year={2025},
  publisher={Oxford University Press}
}
```
