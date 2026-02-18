module BarcodeModels

using ComponentArrays
using DifferentialEquations
using JumpProcesses
using Distributions
using Suppressor
using Random
using DataFrames
using Parameters: @unpack
using StatsBase

include("types/Parameters.jl")
include("types/State.jl")

include("models/abstract.jl")
include("helpers/common_helpers.jl")
include("helpers/ode_helpers.jl")
include("helpers/abm_helpers.jl")
include("models/res_pop.jl")
include("models/res_pop_abm.jl")

include("simulation/simulate_common.jl")
include("simulation/simulate_hybrid.jl")
include("simulation/simulate_abm.jl")
include("simulation/simulate.jl")
include("simulation/noise.jl")

export
    DRUG_EFFECTS,
    ModelParams, ABMParams, ExperimentParams,
    ResPop, ResPop_ABM,
    simulate_experiment

end
