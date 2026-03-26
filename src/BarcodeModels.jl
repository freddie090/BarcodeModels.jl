# BarcodeModels Software Version 0
# Copyright (c) Cancer Research Technology Limited and
# The Institute of Cancer Research. All rights reserved.
# Licensed under the BarcodeModels Academic Use Licence (see LICENSE.txt)

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
include("models/res_dmg.jl")
include("models/res_pop_abm.jl")
include("models/res_dmg_abm.jl")

include("simulation/simulate_common.jl")
include("simulation/simulate_hybrid.jl")
include("simulation/simulate_abm.jl")
include("simulation/simulate.jl")
include("simulation/noise.jl")
include("plotting/simulation_plots.jl")

export
    RESPOP_DRUG_EFFECTS, RESDMG_DRUG_EFFECTS,
    ResPopParams, ResDmgParams,
    ResPopState, ResDmgState,
    ABMParams, ExperimentParams,
    ResPop, ResDmg, ResPop_ABM, ResDmg_ABM,
    simulate_experiment,
    model_meas_noise,
    plot_simulation_outputs

end
