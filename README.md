# BarcodeModels.jl

BarcodeModels.jl is a Julia package for simulating barcode-tracked cell population dynamics under growth, phenotypic switching, and drug treatment using a hybrid ODE + stochastic jump model. It provides a focused API for running long-term barcoded resistance evolution experiments and returning population-level outputs.

## What it does
- Models sensitive/resistant/escape phenotypes with logistic growth and phenotype switching.
- Supports drug effects on birth, death, or combined modes.
- Uses hybrid deterministic/stochastic dynamics to efficiently handle large and small populations.
- Provides experiment-level simulation utilities for expansion and treatment workflows.

## Current package structure
- `src/BarcodeModels.jl`: Package entry point and exports.
- `src/types/Parameters.jl`: Parameter types and validation.
- `src/types/State.jl`: State definitions and helpers.
- `src/models/abstract.jl`: Abstract model type.
- `src/helpers/common_helpers.jl`: Shared helper functions.
- `src/helpers/ode_helpers.jl`: Hybrid/ODE helper functions.
- `src/helpers/abm_helpers.jl`: ABM helper functions.
- `src/models/res_pop.jl`: Hybrid resistance population model core simulation logic.
- `src/models/res_pop_abm.jl`: ABM resistance population model core simulation logic.
- `src/simulation/simulate.jl`: Shared experiment orchestration drivers and API dispatch.
- `test/runtests.jl`: Minimal integration test.

## Quick start
```julia
using BarcodeModels

params = ModelParams(
    b = 1.0,
    d = 0.1,
    rho = 0.0,
    mu = 0.0,
    sig = 0.0,
    del = 0.0,
    al = 0.0,
    Dc = 0.0,
    k = 0.0,
    psi = 0.0,
    drug_effect = :d
)

model = ResPop(params)

exp = ExperimentParams(
    n0 = 10,
    t_exp = 6.0,
    tmax = 10.0,
    t_Pass = -1.0,
    Nseed = 10,
    Nmax = 100,
    Cc = 100,
    treat_ons = Float64[],
    treat_offs = Float64[],
    t_keep = [10.0],
    Nswitch = 100
)

result = simulate_experiment(model, exp)
```

## Key API
- `ModelParams`: Core model parameters (including `rho`, the initial resistant fraction).
- `ResPop`: Resistance population model that validates parameters on construction.
- `ResPop_ABM`: Agent-based resistance population model.
- `ExperimentParams`: Experiment-level configuration (times, passage schedule, sampling).
- `simulate_experiment(model, exp; kwargs...)`: Runs the experiment using the model-class-specific driver.

## Tests
```bash
julia --project -e "using Pkg; Pkg.test()"
```

## Notes
- This package is under active development; API details may change.
- Documentation is minimal at this stage; see source for implementation details.
