module PEtabTraining

using ArgCheck: @argcheck
import Catalyst
import CSV
using DataFrames: DataFrame, nrow
import PEtab: PEtab, PEtabModel, PEtabODEProblem
using SimpleUnPack: @unpack

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "splitting.jl"))
include(joinpath(@__DIR__, "curriculum.jl"))
include(joinpath(@__DIR__, "multiple_shooting.jl"))
include(joinpath(@__DIR__, "petab_odeproblem.jl"))
include(joinpath(@__DIR__, "util.jl"))

export PEtabCurriculumProblem, PEtabMultipleShootingProblem, SplitUniform, SplitCustom

end
