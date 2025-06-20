module PEtabTraining

using ArgCheck: @argcheck
import Catalyst
import CSV
using DataFrames: DataFrame
import PEtab: PEtab, PEtabODEProblem, PEtabModel, get_obs_sd_parameter
using ModelingToolkit
import RuntimeGeneratedFunctions: RuntimeGeneratedFunctions, @RuntimeGeneratedFunction
using SBMLImporter

RuntimeGeneratedFunctions.init(@__MODULE__)

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "splitting.jl"))
include(joinpath(@__DIR__, "curriculum.jl"))
include(joinpath(@__DIR__, "petab_model.jl"))
include(joinpath(@__DIR__, "petab_odeproblem.jl"))

export PEtabCurriculumProblem, SplitUniform, SplitCustom

end
