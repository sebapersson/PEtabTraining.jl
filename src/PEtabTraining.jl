module PEtabTraining

using ArgCheck: @argcheck
import Catalyst
import ComponentArrays
import CSV
import DataFrames: DataFrames, DataFrame, nrow
import PEtab: PEtab, PEtabModel, PEtabODEProblem
using SimpleUnPack: @unpack

include(joinpath(@__DIR__, "structs.jl"))

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "splitting.jl"))
include(joinpath(@__DIR__, "curriculum.jl"))
include(joinpath(@__DIR__, "curriculum_ms_combined.jl"))
include(joinpath(@__DIR__, "multiple_shooting.jl"))
include(joinpath(@__DIR__, "petab_odeproblem.jl"))
include(joinpath(@__DIR__, "util.jl"))

export PEtabClProblem, PEtabMSProblem, PEtabCLMSProblem, SplitTime, SplitData

end
