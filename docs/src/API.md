```@meta
CollapsedDocStrings=true
```

# [API](@id API)

PEtabTraining implements three training strategies that can be applied directly to a
`PEtabODEProblem`:

```@docs; canonical=true
PEtabClProblem
PEtabMsProblem
PEtabClMsProblem
```

All above training strategies require splitting the original simulation in windows/stages:

```@docs; canonical=true
SplitTime
SplitData
```

Utility functions are provided for setting window initial values and the continuity window
penalty for multiple-shooting-based methods.

```@docs; canonical=true
set_ms_window_penalty!
set_u0_ms_windows!
MsInitConstant
MsInitFirst
MsInitSimulate
```

For `PEtabClMsProblem`, there is also a function for mapping parameters between curriculum
stages:

```@docs; canonical=true
map_x_stage
```
