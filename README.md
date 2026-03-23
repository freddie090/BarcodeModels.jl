# BarcodeModels.jl

BarcodeModels.jl is a Julia package for simulating barcoded cell population dynamics under growth, phenotypic switching and drug treatment. The models can reproduce key features of lineage tracing experiments in pre-clinical cancer cell systems, including cell barcoding, a mutual expansion period, sampling into replicate flasks and passaging bottlenecks.

## Model Classes
BarcodeModels is organised around two model **classes**:

- **Hybrid**: deterministic ODE dynamics coupled with stochastic jump dynamics.
- **Agent-Based**: single-cell stochastic simulations.

Model class selection is determined by the model object you pass into `simulate_experiment`.

1. Define biological parameters with model-specific types (e.g. `ResPopParams`).
2. Instantiate a model using either the Hybrid or Agent-Based class.
3. Define the experimental design with `ExperimentParams`.
4. Call `simulate_experiment(model, exp)`.

Dispatch is class-based:
- Passing a Hybrid-class model routes to the hybrid simulation pipeline.
- Passing an Agent-Based-class model routes to the ABM simulation pipeline.


## Current Models

The package contains different models that encompass different hypotheses for how resistance can evolve.

### Shared Model Features

All models in `BarcodeModels.jl` share a common structure describing phenotype-specific population growth under treatment.

#### 1. Drug uptake and decay  
Drug exposure is modelled through a time-dependent effective concentration $D_c(t)$, governed by a simple uptake–decay process:
- `Dc` sets the maximum effective drug concentration ($D_c$)
- `k` controls the rate of drug uptake and/or decay ($\kappa$) 

Treatment is applied through scheduled on/off windows.

#### 2. Phenotype-specific birth, death and transition events  
Cells exist in discrete phenotypic states with each phenotype $X$ defined by:
- Birth rate $b_X$  
- Death rate $d_X$

Population dynamics arise from:
- Birth events (cell division)  
- Death events  
- Phenotype transition events between states  


#### 3. Intrinsic and treatment-dependent transitions  
Phenotype transitions can occur via two distinct mechanisms:
- **Intrinsic transitions**: occur independently of treatment and reflect baseline plasticity  
- **Treatment-dependent transitions**: occur as a function of drug exposure,  modulated by $D_c(t)$  

This allows the models to capture both pre-existing and drug-induced resistance mechanisms.

#### 4. Logistic population growth  
Population expansion is regulated by a carrying capacity $C_c$, implemented through logistic growth.

### ResPop
`ResPop` models up to three possible phenotypes (`S`, `R`, `E`) and can be run using either model class:
- Hybrid class: `ResPop(params)`
- Agent-Based class: `ResPop_ABM(params)`

#### ResPop features:
- Phenotypes: sensitive (`S`), resistant (`R`), and escape (`E`).
- Phenotype transitions: `S→R`, `R→S`, `R→E`.
- Transition behaviour:
    - `S→R` and `R→S` are intrinsic transitions that are independent of treatment.
    - `R→E` occurs as a function of effective treatment concentration to capture treatment-induced escape transitions.
- Impact of treatment on phenotype birth/death rates:
    - All phenotypes are affected by treatment through the selected drug-effect mode (`:d`, `:b`, `:c`).
        - `:d`: drug increases the death rate
        - `:b`: drug decreases the birth rate
        - `:c`: drug decreases the birth rate and then increases the death rate
    - Phenotypes `R` and `E` experience reduced effect of the drug
        - e.g. for `:d` drug effect: $d_R = d + (D_C(t) \cdot (1 - \psi))$

### ResDmg
`ResDmg` models up to four possible phenotypes (`S`, `D`, `R`, `E`) and can be run using either model class:
- Hybrid class: `ResDmg(params)`
- Agent-Based class: `ResDmg_ABM(params)`

#### ResDmg features:
- States: sensitive (`S`), damaged (`D`), resistant (`R`), and escape (`E`).
- Phenotype transitions: `S→R`, `R→S`, `R→E`, `S→D`, `R→D`, `D→S`.
- Transition behaviour:
    - `S→R` and `R→S` are intrinsic transitions that are independent of treatment.
    - `S→D` and `R→D` capture treatment-induced transitions into the damaged states.
    - `D→S` (repair) is available once cells are damaged and represents intrinsic recovery dynamics.
    - `R→E` occurs as a function of effective treatment concentration to capture treatment-induced escape transitions.
- Impact of treatment on phenotype birth/death rates:
    - Treatment can increase damage-state occupancy through `S→D` and `R→D` routes.
    - Damaged cells do not divide and experience elevated death rates.
    - The resistant phenotype (`R`) experiences a lower probability of drug-induced damage via `psi`.




## Experimental design

`BarcodeModels.jl` uses a shared **experimental design layer**, so you can run different model implementations with the same `ExperimentParams` object.

This allows different models to be compared under the same fixed experimental parameters.

Features of the experimental design include:
- Shared `ExperimentParams` across all model families, enabling direct model-to-model comparisons.
- Support for scalar or staged expansion (`t_exp`) with matching bottleneck seeding rules (`Nseed`).
- Scheduled passaging by time (`t_Pass`) and/or cap-triggered passaging (`Nmax`).
- Replicate experiments (`n_rep`) under identical fixed design settings.
- Optional treatment windows (`treat_ons`, `treat_offs`) and compact observation outputs (`t_keep`).
- Optional assay branches (`run_IC`, `run_colony`) where enabled in simulation pipelines.


## Quick start

### 1) Download or clone the code

```bash
git clone https://github.com/freddie090/BarcodeModels.jl.git
cd BarcodeModels.jl
```

### 2) Activate the project environment and install dependencies

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### 3) Load the package

```julia
using BarcodeModels
```

### 4) Run a minimal experiment

```julia
params = ResPopParams(
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
    t_Pass = Float64[],
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

# ResDmg variant
params_dmg = ResDmgParams(
    b = 0.893,
    d = 0.200,
    rho = 1e-01,
    mu = 1e-01,
    sig = 0.0,
    del = 0.0,
    ome = 10.0,
    zet_S = 0.01,
    zet_R = 0.01,
    Dc = 1.4,
    k = 0.5,
    psi = 0.0,
    drug_effect = :c
)

model_dmg_hybrid = ResDmg(params_dmg)
model_dmg_abm = ResDmg_ABM(params_dmg)

result_dmg_hybrid = simulate_experiment(model_dmg_hybrid, exp)
result_dmg_abm = simulate_experiment(model_dmg_abm, exp)
```

## Parameter reference

### `ResPopParams`
Core parameters for the ResPop model family (hybrid and ABM).

| Parameter | Type (default) | Meaning | Constraints |
|---|---|---|---|
| `b` | `Float64` (required) | Baseline (Sensitive) per-cell birth rate ($b$). | Must be valid numeric input. |
| `d` | `Float64` (required) | Baseline (Sensitive) per-cell death rate ($d$). | Must be valid numeric input. |
| `rho` | `Float64` (`0.0`) | Initial resistant fraction at initial barcoding `(t=0)` ($\rho$). | `0 ≤ rho ≤ 1`. |
| `mu` | `Float64` (required) | `S→R` switching probability per division ($\mu$). | `0 ≤ mu ≤ 1`. |
| `sig` | `Float64` (required) | `R→S` switching probability per division ($\sigma$). | `0 ≤ sig ≤ 1`. |
| `del` | `Float64` (required) | Resistant fitness-cost term ($\delta$). | `0 ≤ del ≤ 1`. |
| `al` | `Float64` (required) | `R→E` switching probability per division ($\alpha$). | `0 ≤ al ≤ 1`, and `al + sig ≤ 1`. |
| `Dc` | `Float64` (required) | Drug effect magnitude scale ($D_c$). | If `drug_effect == :b`, must satisfy additional bounds (see below). |
| `k` | `Float64` (required) | Drug concentration accumulation/decay rate ($\kappa$). | Must be valid numeric input. |
| `psi` | `Float64` (required) | Drug attenuation factor for resistant/escape phenotypes ($\psi$). | `0 ≤ psi ≤ 1`. |
| `drug_effect` | `Symbol` (`:d`) | Drug action mode: death-only (`:d`), birth-only (`:b`), combined (`:c`). | Must be one of `:d`, `:b`, `:c`. |

### `ResDmgParams`
Core parameters for the ResDmg model family (hybrid and ABM). ResDmg uses `S, DS, DR, R` states and does not include an escape (`E`) transition.

| Parameter | Type (default) | Meaning | Constraints |
|---|---|---|---|
| `b, d, rho, mu, sig, del, Dc, k, psi, drug_effect` | as in `ResPopParams` (except `al`) | Same meanings as ResPop parameters where shared ($b, d, \rho, \mu, \sigma, \delta, D_c, \kappa, \psi$). | Shared constraints from `ResPopParams` where applicable. |
| `ome` | `Float64` (required) | Damage transition rate ($\omega$) into damaged states (`S->DS`, `R->DR`, with `R` damage attenuated by `1-psi`). | `ome >= 0`. |
| `zet_S` | `Float64` (required) | Repair transition rate ($\zeta_S$) from `DS->S`. | `zet_S >= 0`. |
| `zet_R` | `Float64` (required) | Repair transition rate ($\zeta_R$) from `DR->R`. | `zet_R >= 0`. |

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

## Key API
- `ResPopParams`: ResPop biological parameters.
- `ResDmgParams`: ResDmg biological parameters.
- `ResPop`: hybrid ODE+jump model.
- `ResDmg`: hybrid ODE+jump model with damaged compartment.
- `ResPop_ABM`: agent-based model.
- `ResDmg_ABM`: agent-based model with damaged compartment.
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
- `src/models/res_dmg.jl`: hybrid ResDmg model core simulation logic.
- `src/models/res_pop_abm.jl`: ABM model core simulation logic.
- `src/models/res_dmg_abm.jl`: ABM ResDmg model core simulation logic.
- `src/simulation/simulate.jl`: experiment orchestration drivers and API dispatch.
- `test/runtests.jl`: integration tests.

## Tests
```bash
julia --project -e "using Pkg; Pkg.test()"
```
