module PEtabTraining

using ArgCheck: @argcheck
import ComponentArrays: ComponentArrays, ComponentArray
import CSV
import DataFrames: DataFrames, DataFrame, nrow
import PEtab: PEtab, PEtabModel, PEtabODEProblem
using SimpleUnPack: @unpack
import StatsBase
using StyledStrings: styled, @styled_str
using Printf: @sprintf

include(joinpath(@__DIR__, "structs.jl"))

include(joinpath(@__DIR__, "cl_ms_combined.jl"))
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "curriculum.jl"))
include(joinpath(@__DIR__, "multiple_shooting.jl"))
include(joinpath(@__DIR__, "petab_odeproblem.jl"))
include(joinpath(@__DIR__, "splitting.jl"))
include(joinpath(@__DIR__, "show.jl"))
include(joinpath(@__DIR__, "util.jl"))

export PEtabClProblem, PEtabMsProblem, PEtabClMsProblem, SplitTime, SplitData,
    MsInitConstant, MsInitFirst, MsInitSimulate, set_u0_ms_windows!, set_ms_window_penalty!,
    map_x_stage, allocate_cl_epochs

end
