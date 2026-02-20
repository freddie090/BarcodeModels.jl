# BarcodeModels.jl

BarcodeModels.jl is a Julia package for simulating barcoded cell population dynamics under growth, phenotypic switching, passage bottlenecks, and drug treatment.

## Model families
BarcodeModels provides two model implementations of the same biological process under the same experiment design interface (`ExperimentParams` + `simulate_experiment`):

- `ResPop` (hybrid ODE + stochastic jumps)
- `ResPop_ABM` (agent-based model)

The hybrid model is more memory efficient and faster than the agent-based model, while both can simulate the same process with the same experimental design settings.

## What it does
- Models sensitive (`S`), resistant (`R`), and escape (`E`) phenotypes with logistic population constraint.
- Supports phenotype switching (`S→R`, `R→S`, `R→E`) and treatment schedules.
- Supports drug effects on death (`:d`), birth (`:b`), or combined (`:c`) dynamics.
- Supports expansion and bottleneck passage workflows with repeated replicates.

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

# Choose either model family
model_hybrid = ResPop(params)
model_abm = ResPop_ABM(params)

exp = ExperimentParams(
    n0 = 10,
    t_exp = 6.0,
    tmax = 10.0,
    t_Pass = Float64[],  # no scheduled passage times
    Nseed = 10,
    Nmax = 100,
    Cc = 100,
    treat_ons = Float64[],
    treat_offs = Float64[],
    t_keep = [10.0],
    Nswitch = 100
)

result_hybrid = simulate_experiment(model_hybrid, exp)
result_abm = simulate_experiment(model_abm, exp)
```

## Parameter reference

### `ModelParams`
Core biological/process parameters shared by hybrid and ABM models.

| Parameter | Type (default) | Meaning | Constraints |
|---|---|---|---|
| `b` | `Float64` (required) | Baseline per-cell birth rate. | Must be valid numeric input. |
| `d` | `Float64` (required) | Baseline per-cell death rate. | Must be valid numeric input. |
| `rho` | `Float64` (`0.0`) | Initial resistant fraction at seeding. | `0 ≤ rho ≤ 1`. |
| `mu` | `Float64` (required) | `S→R` switching coefficient. | `0 ≤ mu ≤ 1`. |
| `sig` | `Float64` (required) | `R→S` switching coefficient. | `0 ≤ sig ≤ 1`. |
| `del` | `Float64` (required) | Resistant fitness cost term. | `0 ≤ del ≤ 1`. |
| `al` | `Float64` (required) | `R→E` switching coefficient. | `0 ≤ al ≤ 1`, and `al + sig ≤ 1`. |
| `Dc` | `Float64` (required) | Drug effect magnitude scale. | If `drug_effect == :b`, must satisfy additional bounds (see below). |
| `k` | `Float64` (required) | Drug concentration rise/decay rate control. | Must be valid numeric input. |
| `psi` | `Float64` (required) | Drug attenuation factor for resistant/escape phenotypes. | `0 ≤ psi ≤ 1`. |
| `drug_effect` | `Symbol` (`:d`) | Drug action mode: death-only, birth-only, combined. | Must be one of `:d`, `:b`, `:c`. |

Additional `drug_effect == :b` requirement:
- `Dc <= b`
- if `psi < 1`, `Dc*(1-psi) <= b*(1-del)`

### `ABMParams`
ABM runtime and barcode-library controls.

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
