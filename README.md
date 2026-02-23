# BarcodeModels.jl

BarcodeModels.jl is a Julia package for simulating barcoded cell population dynamics under growth, phenotypic switching and drug treatment. The models can reproduce key features of lineage tracing experiments in pre-clinical cancer cell systems, including cell barcoding, a mutual expansion period, sampling into replicate flasks and passaging bottlenecks.

## Model types
BarcodeModels is organized around two model **types**:

- **Hybrid models**: deterministic ODE dynamics coupled with stochastic jump dynamics.
- **Agent-based models (ABM)**: explicit single-cell stochastic simulations.

Both model types can simulate the same biological process under the same experimental design settings.

## Choosing a Model
Model selection is determined by the model object you pass into `simulate_experiment`.

1. Define shared biological parameters with `ModelParams`.
2. Instantiate a model class belonging to a model type (hybrid or ABM).
3. Define experimental design with `ExperimentParams`.
4. Call `simulate_experiment(model, exp)`.

Dispatch is type-based:
- Passing a hybrid model class routes to the hybrid simulation pipeline.
- Passing an ABM model class routes to the ABM simulation pipeline.

## Current model classes
The current package versions of these model types are:

- `ResPop` (hybrid model class)
- `ResPop_ABM` (agent-based model class)

## ResPop model features
- Models sensitive (`S`), resistant (`R`), and escape (`E`) phenotypes with logistic population constraint.
- Supports phenotype switching (`S→R`, `R→S`, `R→E`) with probabilities coupled to cell division events.
- Uses `Nswitch` for phenotype birth/death ODE↔jump switching and `N_trans_switch` for transition-process ODE↔jump switching.
- Supports drug effects on death (`:d`), birth (`:b`), or combined (`:c`) dynamics.
- Supports expansion and bottleneck passage workflows with experimental replicates.

## Experimental design

`BarcodeModels.jl` uses a shared **experimental design layer**, so you can run different model implementations with the same `ExperimentParams` object.

This allows different models to be compared under the same fixed experimental parameters.


## Quick start
```julia
using BarcodeModels

params = ModelParams(
    b = 0.893,
    d = 0.200,
    rho = 1e-01,
    mu = 1e-01,
    sig = 0.0,
    del = 0.0,
    al = 0.0,
    Dc = 1.4,
    k = 0.5,
    psi = 0.0,
    drug_effect = :c
)

# Choose either model class (from hybrid or ABM model types)
model_hybrid = ResPop(params)
model_abm = ResPop_ABM(params)

exp = ExperimentParams(
    n0 = Int(1e+03),
    t_exp = 6.0,
    tmax = 30.0,
    t_Pass = Float64[],  # no scheduled passage times
    Nseed = Int(1e+03),
    Nmax = Int(5e+04),
    Cc = Int(1e+05),
    treat_ons = Float64[1.0],
    treat_offs = Float64[30.0],
    t_keep = [5.0],
    Nswitch = 500,
    full_sol = true,
    n_rep = 4
)

result_hybrid = simulate_experiment(model_hybrid, exp)
result_abm = simulate_experiment(model_abm, exp)
```

## Parameter reference

### `ModelParams`
Core biological/process parameters shared by hybrid and ABM models.

| Parameter | Type (default) | Meaning | Constraints |
|---|---|---|---|
| `b` | `Float64` (required) | Baseline (Sensitive) per-cell birth rate. | Must be valid numeric input. |
| `d` | `Float64` (required) | Baseline (Sensitive) per-cell death rate. | Must be valid numeric input. |
| `rho` | `Float64` (`0.0`) | Initial resistant fraction at initial barcoding `(t=0)` | `0 ≤ rho ≤ 1`. |
| `mu` | `Float64` (required) | `S→R` switching probability (per cell division). | `0 ≤ mu ≤ 1`. |
| `sig` | `Float64` (required) | `R→S` switching probability (per cell division). | `0 ≤ sig ≤ 1`. |
| `del` | `Float64` (required) | Resistant fitness cost term. | `0 ≤ del ≤ 1`. |
| `al` | `Float64` (required) | `R→E` switching probability (per cell division - currently scaled by effective drug concentration). | `0 ≤ al ≤ 1`, and `al + sig ≤ 1`. |
| `Dc` | `Float64` (required) | Drug effect magnitude scale. | If `drug_effect == :b`, must satisfy additional bounds (see below). |
| `k` | `Float64` (required) | Drug concentration accumulation/decay rate. | Must be valid numeric input. |
| `psi` | `Float64` (required) | Drug attenuation factor for resistant/escape phenotypes. | `0 ≤ psi ≤ 1`. |
| `drug_effect` | `Symbol` (`:d`) | Drug action mode: death-only (`:d`), birth-only (`:b`), combined (`:c`). | Must be one of `:d`, `:b`, `:c`. |

\
Additional `drug_effect == :b` requirement:
- `Dc <= b`
- if `psi < 1`, `Dc*(1-psi) <= b*(1-del)`
\

### `ABMParams`
Parameters specific to the agent-based model implementation.

| Parameter | Type (default) | Meaning | Constraints / notes |
|---|---|---|---|
| `Nbuff` | `Int64` (`100000`) | Cell-buffer size (live + dead slots). | Increase for large populations; affects memory. |
| `t_frac` | `Float64` (`0.005`) | Fraction of `tmax` used for ABM recording interval. | Must be in `[0, 1]` at runtime. |
| `dt_save_at` | `Float64` (`0.1`) | Drug concentration interpolation/save cadence for ABM. | Positive recommended. |
| `skew_lib` | `Bool` (`false`) | Enable skewed barcode library initialization. | ABM-only feature. |
| `bc_unif` | `Float64` (`0.0`) | Skew/uniformity control for barcodes. | Used when `skew_lib=true`. |
| `Nbc` | `Int64` (`0`) | Number of barcodes in skewed library mode. | Used when `skew_lib=true`. |
| `sub_sample_cells` | `Bool` (`false`) | Enable per-passage subsampling outputs. | ABM-only output option. |
| `K` | `Int64` (`0`) | Subsample size when `sub_sample_cells=true`. | Must be consistent with live-cell count. |

### `ExperimentParams`
Experiment design shared by both model families.

| Parameter | Type (default) | Meaning | Constraints / notes |
|---|---|---|---|
| `n0` | `Int64` (required) | Initial seeding population. | `> 0`. |
| `t_exp` | `Float64` or `Vector{Float64}` (required) | Expansion phase duration(s). | Scalar or vector; if vector, must align with `Nseed` vector usage. |
| `tmax` | `Float64` (required) | Final experiment time horizon. | `> 0`. |
| `t_Pass` | `Float64` or `Vector{Float64}` (required) | Scheduled passage times. | Strictly `> 0`; use `Float64[]` for no passage events. |
| `Nseed` | `Int64` or `Vector{Int64}` (required) | Passage bottleneck sample size(s). | Positive; vector form must align with expansion design. |
| `Nmax` | `Int64` (required) | Population cap trigger for growth/passaging logic. | `> 0` and `Nmax <= Cc`. |
| `Cc` | `Int64` (required) | Carrying-capacity scale in logistic terms. | `> 0`. |
| `treat_ons` | `Vector{Float64}` (required) | Treatment on-times. | Vector of numeric times. |
| `treat_offs` | `Vector{Float64}` (required) | Treatment off-times. | Vector of numeric times. |
| `t_keep` | `Vector{Float64}` (required) | Observation times to retain in compact outputs. | Vector of numeric times. |
| `Nswitch` | `Int64` (required) | Hybrid model ODE/jump switching threshold. | `> 0`. |
| `N_trans_switch` | `Float64` (`1000.0`) | Hybrid transition-activity threshold for deterministic vs jump switching (`S↔R`, `R→E`). | `> 0`. |
| `save_at` | `Float64` (`0.5`) | Hybrid solver save interval. | `> 0`. |
| `n_rep` | `Int64` (`4`) | Number of replicate experiments. | `> 0`. |
| `drug_treatment` | `Bool` (`true`) | Whether treatment is applied in main experiment runs. | Boolean. |
| `full_sol` | `Bool` (`false`) | Return full per-time simulated trajectory table. | Boolean. |
| `run_IC` | `Bool` (`false`) | Run initial-condition assay branch. | Boolean. |
| `IC_n0` | `Int64` (`1000`) | Initial-condition assay seed size. | `> 0`. |
| `IC_tmax` | `Float64` (`4.0`) | Initial-condition assay end time. | `> 0`. |
| `IC_treat_on` | `Float64` (`1.0`) | Initial-condition treatment start time. | Numeric. |
| `run_colony` | `Bool` (`false`) | Run colony-formation assay branch. | Boolean. |
| `nCol` | `Int64` (`1000`) | Number of colony simulations. | `> 0`. |
| `tCol` | `Float64` (`12.0`) | Colony assay end time. | `> 0`. |
| `ColNmax` | `Int64` (`50`) | Colony endpoint threshold for "successful" growth. | `> 0`. |

## Key API
- `ModelParams`: shared biological parameters.
- `ResPop`: hybrid ODE+jump model.
- `ResPop_ABM`: agent-based model.
- `ABMParams`: ABM-specific configuration.
- `ExperimentParams`: experiment-level design parameters.
- `simulate_experiment(model, exp; kwargs...)`: run simulation for either model family.

## Current package structure
- `src/BarcodeModels.jl`: package entry point and exports.
- `src/types/Parameters.jl`: parameter types and validation.
- `src/types/State.jl`: state definitions and helpers.
- `src/models/abstract.jl`: abstract model type.
- `src/helpers/common_helpers.jl`: shared helper functions.
- `src/helpers/ode_helpers.jl`: hybrid/ODE helper functions.
- `src/helpers/abm_helpers.jl`: ABM helper functions.
- `src/models/res_pop.jl`: hybrid model core simulation logic.
- `src/models/res_pop_abm.jl`: ABM model core simulation logic.
- `src/simulation/simulate.jl`: experiment orchestration drivers and API dispatch.
- `test/runtests.jl`: integration tests.

## Tests
```bash
julia --project -e "using Pkg; Pkg.test()"
```
