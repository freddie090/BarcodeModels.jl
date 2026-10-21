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

`ResPopInVivo` is available in both classes, extending ResPop with in vivo-specific engraftment parameters and implantation/passaging behaviour. See [InVivo functionality](#invivo-functionality) for details.


## Current Models

The package contains different models that encompass different hypotheses for how resistance can evolve.

### Shared Model Features

All models in `BarcodeModels.jl` share a common structure describing phenotype-specific population growth under treatment.

Treatment is applied through scheduled on/off windows.

#### 1. Phenotype-specific birth, death and transition events  
Cells exist in discrete phenotypic states with each phenotype $X$ defined by:
- Birth rate $b_X$  
- Death rate $d_X$

Population dynamics arise from:
- Birth events (cell division)  
- Death events  
- Phenotype transition events between states  

#### 2. Intrinsic and treatment-dependent transitions  
Phenotype transitions can occur via two distinct mechanisms:
- **Intrinsic transitions**: occur independently of treatment
- **Treatment-dependent transitions**: occur as a function of drug exposure,  modulated by the current effective concentration ($D_c(t)$)

This allows the models to capture both pre-existing and drug-induced resistance mechanisms.

#### 3. Logistic population growth  
Population growth is logistic and is regulated by a carrying capacity $C_c$.

#### 4. Drug uptake and decay  
Drug exposure is modelled through a time-dependent effective concentration $D_c(t)$, governed by a simple uptake–decay process:
- `Dc` sets the maximum effective drug concentration ($D_c$)
- `k` controls the rate of drug uptake and/or decay ($\kappa$) 


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
    - Phenotypes `R` and `E` experience reduced effect of the drug mediated via `psi` (where `psi=0` means no attenuation and `psi=1` means complete attenuation).
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




## Simple simulations (single population)

`BarcodeModels.jl` also provides a simplified **single-population simulation layer** through `simulate_simple(model, sim)`.

This method is intended for quick model implementation without the full replicate/passaging experiment workflow. This version is useful for testing and debugging and for generating single-run outputs for a given model and parameter set.

Features of the simple simulation layer include:
- Shared dispatch across model families (`ResPop`, `ResDmg`, InVivo variants, and EvBC ABM variants).
- Input via `SimpleSimParams` for one-population runs with optional treatment windows.
- Compact outputs for trajectory and barcode summaries (`sol_df`; plus `lin_df` for ABM classes).
- EvBC support returns lineage outputs (`lineage_df`) with ground-truth relatedness metadata.
- No experiment-level expansion/passaging/replicate orchestration.

### Outputs

- `sol_df`: time-series population output table for the simulated run.
- `lin_df`: barcode abundance table for ABM model classes.
- `lineage_df`: EvBC-only ground-truth lineage table with parent-child relationships and lineage metadata.

InVivo simple simulations return EG-stratified population counts and support optional ABM phenotype-by-barcode outputs. Engraftment, control replicates and pool splitting belong to the full experiment workflow.

Runnable examples are provided in `Quick start`.


## Experimental design

`BarcodeModels.jl` uses a shared **experimental design layer**, so you can run different model implementations with the same `ExperimentParams` object.

This allows different models to be compared under the same fixed experimental parameters.

Features of the experimental design include:
- Shared `ExperimentParams` across all model families, enabling direct model-to-model comparisons.
- Support for scalar or staged expansion (`t_exp`) with matching bottleneck seeding rules (`Nseed`).
- Scheduled passaging by time (`t_Pass`) and/or cap-triggered passaging (`Nmax`).
- Replicate experiments (`n_rep`) under identical fixed design settings.
- Optional treatment windows (`treat_ons`, `treat_offs`) and compact observation outputs (`t_keep`).


### Outputs

- `sol_df`: full experiment trajectory table across time and replicate outputs.
- `lin_df`: barcode abundance summary table across replicate/passage outputs.
- `t`: summary vector of sampled times (returned unless `just_lin=true`).
- `u`: summary vector of population sizes aligned to `t` (returned unless `just_lin=true`).
- `sub_lin_df`: optional subsampled barcode table when `sub_sample_cells=true`.
- `lineage_df`: EvBC-only ground-truth lineage table with node-level lineage metadata.

InVivo experiments also return engraftment counts and condition labels, with endpoint summaries for each condition/replicate. See [InVivo outputs](#invivo-outputs-and-how-to-access-them) for their output conventions.


## Quick start

### 1) Install from GitHub (Julia package manager)

```julia
import Pkg
Pkg.add(url = "https://github.com/freddie090/BarcodeModels.jl.git")
```

### 2) Update to the latest package version

```julia
import Pkg
Pkg.update("BarcodeModels")
```

### 3) Load the package

```julia
using BarcodeModels
```

### 4) Run a minimal simple simulation (single population)

```julia
# Define the parameters describing the biological model and treatment effects
# (biological model defined by the type of params object used)
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

# Define the parameters describing the experimental variables
sim = SimpleSimParams(
    n0 = Int(1e+03),
    tmax = 30.0,
    Nmax = Int(5e+04),
    Cc = Int(1e+05),
    treat_ons = Float64[1.0],
    treat_offs = Float64[30.0],
    Nswitch = 500,
    Aswitch = 1000.0,
    save_at = 0.5,
    drug_treatment = true
)

# Simulate simple: Hybrid
simple_hybrid = simulate_simple(ResPop(params), sim)


```

### 5) Run a minimal full experiment

```julia
# Define the parameters describing the biological model and treatment effects
# (biological model defined by the type of params object used)
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
abm = ABMParams(Nbuff = 100000, t_frac = 0.005, dt_save_at = 0.1)
model_abm = ResPop_ABM(params; abm = abm)

# Define the parameters describing the experimental variables
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

# Simulate experiment: Hybrid
result_hybrid = simulate_experiment(model_hybrid, exp)
# Simulate experiment: ABM
result_abm = simulate_experiment(model_abm, exp)

# ResDmg model variant
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

# Choose either model class (from hybrid or ABM model types)
model_dmg_hybrid = ResDmg(params_dmg)
model_dmg_abm = ResDmg_ABM(params_dmg)

# Simulate experiment: Hybrid
result_dmg_hybrid = simulate_experiment(model_dmg_hybrid, exp)
# Simulate experiment: ABM
result_dmg_abm = simulate_experiment(model_dmg_abm, exp)
# Note that the same `exp` object can be used for both model classes
```

### Phenotypic Transition Modes

All hybrid, ABM, and EvBC model constructors support both modes of phenotypic switching:
- **Birth-coupled:** *transitions occur only at division events, with probabilities per division. One daughter cell inherits the parent's phenotype. **(Default)***
- **Independent:** *transitions occur independently of division, with rates per day.*

Usage example:

```julia
coupled = ResPopParams(b=0.8, d=0.2, mu=0.01, sig=0.02, del=0.1,
    al=0.03, Dc=0.4, k=1.0, psi=0.8, transition_coupling=:birth_coupled)
independent = ResPopParams(b=0.8, d=0.2, mu=0.01, sig=0.02, del=0.1,
    al=0.03, Dc=0.4, k=1.0, psi=0.8, transition_coupling=:independent_rate)
```



## Parameter reference

All ResPop, ResDmg and ResPopInVivo model constructors require `b > d` for positive baseline growth. Effective growth under treatment can still be negative.

### `ResPopParams`
Core parameters for the ResPop model family (hybrid and ABM).

| Parameter | Type (default) | Meaning | Constraints |
|---|---|---|---|
| `b` | `Float64` (required) | Baseline (Sensitive) per-cell birth rate ($b$). | Must be valid numeric input; `b > d`. |
| `d` | `Float64` (required) | Baseline (Sensitive) per-cell death rate ($d$). | Must be valid numeric input; `d < b`. |
| `rho` | `Float64` (`0.0`) | Initial resistant fraction at initial barcoding `(t=0)` ($\rho$). | `0 ≤ rho ≤ 1`. |
| `mu` | `Float64` (required) | S division outcome probability in coupled mode; S→R rate/day in independent mode ($\mu$). | Birth coupled: `0 ≤ mu ≤ 1`; independent: `0 ≤ mu`. |
| `sig` | `Float64` (required) | R division outcome probability in coupled mode; R→S rate/day in independent mode ($\sigma$). | Birth coupled: `0 ≤ sig ≤ 1`; independent: `0 ≤ sig`. |
| `del` | `Float64` (required) | Resistant fitness-cost term ($\delta$). | `0 ≤ del ≤ 1`. |
| `al` | `Float64` (required) | R→E daughter outcome factor in coupled mode; R→E rate/day in independent mode ($\alpha$). Both multiply current treatment activation γ. | Birth coupled: `0 ≤ al ≤ 1`, `al + sig ≤ 1`; independent: `0 ≤ al`. |
| `Dc` | `Float64` (required) | Drug effect magnitude scale ($D_c$). | If `drug_effect == :b`, must satisfy additional bounds (see below). |
| `k` | `Float64` (required) | Drug concentration accumulation/decay rate ($\kappa$). | Must be valid numeric input. |
| `psi` | `Float64` (required) | Drug attenuation factor for resistant/escape phenotypes ($\psi$). | `0 ≤ psi ≤ 1`. |
| `drug_effect` | `Symbol` (`:d`) | Drug action mode: death-only (`:d`), birth-only (`:b`), combined (`:c`). | Must be one of `:d`, `:b`, `:c`. |
| `transition_coupling` | `Symbol` (`:birth_coupled`) | Division outcomes or independent conversion of existing cells. | `:birth_coupled` or `:independent_rate`; exact string equivalents also accepted. |

`R_real` controls how the resistant fitness cost `del` is realised in **all Hybrid, ABM, and EvBC models**, including ResPopInVivo. The default is `"l"` in every class. Set it when calling either simulation API:

```julia
result = simulate_experiment(model, exp; R_real = "l")
simple = simulate_simple(model, sim; R_real = "d")
```

| `R_real` | Resistant birth rate | Resistant death rate | Meaning |
|---|---|---|---|
| `"l"` (default) | `b*(1-del)` | `d*(1-del)` | Reduce both birth and death. |
| `"b"` | `b-(b-d)*del` | `d` | Reduce birth only. |
| `"d"` | `b` | `d+(b-d)*del` | Increase death only. |

These are undamaged R rates before treatment, density scaling, and phenotype transitions. All three give the same baseline birth-minus-death rate, `(b-d)*(1-del)`, but different turnover and, with birth-coupled switching, different daughter-production rates. ResPop E uses baseline birth/death rates; ResDmg damaged states do not divide and have death rate `d + Dc` before density scaling. The selected option applies throughout expansion, replicate growth, and passaging.

### `ResPopInVivoParams`

Core parameters for the ResPopInVivo model family (hybrid and ABM). ResPopInVivo uses the ResPop biological parameters plus a fixed engraftment phenotype (`EG0` or `EG1`).

| Parameter | Type (default) | Meaning | Constraints |
|---|---|---|---|
| `b, d, rho, mu, sig, del, al, Dc, k, psi, drug_effect, transition_coupling` | as in `ResPopParams` | Same growth, treatment and phenotype-transition meanings as ResPop. | Shared constraints from `ResPopParams`; `transition_coupling = :birth_coupled` by default. |
| `fEG1` | `Float64` (required) | Initial fraction with the EG1 phenotype, assigned independently of resistance. | Finite and `0 ≤ fEG1 ≤ 1`; a cell fraction in normal expansion, or an observed-barcode fraction with `split_after_barcoding=true`. |
| `pEG` | `Float64` (required) | Probability that an EG1 cell survives implantation or passage engraftment. | Finite and `0 ≤ pEG ≤ 1`. |
| `sEG` | `Float64` (required) | Relative engraftment disadvantage of EG0 cells, whose survival probability is `pEG*(1-sEG)`. | Finite and `0 ≤ sEG ≤ 1`. |

With `split_after_barcoding=true`, `rho` also refers to a fraction of unique observed barcodes rather than cells. The engraftment phenotype is inherited at birth and remains fixed during S/R/E transitions.

### `ResDmgParams`
Core parameters for the ResDmg model family (hybrid and ABM). ResDmg uses `S, DS, DR, R` states and does not include an escape (`E`) transition.

| Parameter | Type (default) | Meaning | Constraints |
|---|---|---|---|
| `b, d, rho, mu, sig, del, Dc, k, psi, drug_effect` | as in `ResPopParams` (except `al`) | Same meanings as ResPop parameters where shared ($b, d, \rho, \mu, \sigma, \delta, D_c, \kappa, \psi$). | Shared constraints from `ResPopParams` where applicable. |
| `ome` | `Float64` (required) | Damage transition rate ($\omega$) into damaged states (`S->DS`, `R->DR`, with `R` damage attenuated by `1-psi`). | `ome >= 0`. |
| `zet_S` | `Float64` (required) | Repair transition rate ($\zeta_S$) from `DS->S`. | `zet_S >= 0` and `zet_R >= zet_S`. |
| `zet_R` | `Float64` (required) | Repair transition rate ($\zeta_R$) from `DR->R`. | `zet_R >= 0` and `zet_R >= zet_S`. |
| `transition_coupling` | `Symbol` (`:birth_coupled`) | As for ResPop, applying to undamaged S↔R only. | `:birth_coupled` or `:independent_rate`; exact string equivalents also accepted. |

\
Additional `drug_effect == :b` requirement:
- `Dc <= b`
- if `psi < 1`, `Dc*(1-psi)` must not exceed the resistant birth rate in the `R_real` table above. This option-dependent bound is checked when setting up the simulation; `Dc <= b` is also checked by model constructors.
\



### `ABMParams`
Parameters specific to the agent-based model implementation.

Barcode initialization supports three modes: default sequential IDs, generated skewed libraries (`skew_lib=true`), or caller-supplied library probabilities (`use_lib_probs=true`). The two library modes are mutually exclusive, and `bc_probs` must already be a valid probability vector with length `Nbc`.

| Parameter | Type (default) | Meaning | Constraints / notes |
|---|---|---|---|
| `Nbuff` | `Int64` (`100000`) | Cell-buffer size (live + dead slots). | Increase for large populations; affects memory. |
| `t_frac` | `Float64` (`0.005`) | Fraction of `tmax` used for ABM recording interval. | Must be in `[0, 1]` at runtime. |
| `dt_save_at` | `Float64` (`0.1`) | Drug concentration interpolation/save cadence for ABM. | Positive recommended. |
| `skew_lib` | `Bool` (`false`) | Enable skewed barcode library initialization. | ABM-only feature. |
| `use_lib_probs` | `Bool` (`false`) | Enable caller-supplied barcode library probabilities. | Mutually exclusive with `skew_lib`. |
| `bc_unif` | `Float64` (`0.0`) | Skew/uniformity control for barcodes. | Used when `skew_lib=true`. |
| `Nbc` | `Int64` (`0`) | Number of barcodes in the library. | Used when `skew_lib=true` or `use_lib_probs=true`. |
| `bc_probs` | `Vector{Float64}` (`Float64[]`) | Library barcode probabilities used for direct barcode sampling. | Used when `use_lib_probs=true`; must have length `Nbc` and sum to 1. |
| `sub_sample_cells` | `Bool` (`false`) | Enable per-passage subsampling outputs. | ABM-only output option. |
| `K` | `Int64` (`0`) | Subsample size when `sub_sample_cells=true`. | Must be consistent with live-cell count. |
| `split_after_barcoding` | `Bool` (`false`) | Split the barcoded starting pool directly into InVivo replicates, without expansion. | InVivo ABM experiments only; requires `use_lib_probs=true` and valid `Nbc`/`bc_probs`. |
| `full_pheno_bc` | `Bool` (`false`) | Return S/R/E counts for each barcode in `pheno_bc_df`. | InVivo ABM only; counts combine EG0 and EG1. |

### `ExperimentParams`
Experiment design shared by the model families. The `inc_control` and `inc_pot` settings apply only to InVivo experiments.

| Parameter | Type (default) | Meaning | Constraints / notes |
|---|---|---|---|
| `n0` | `Int64` (required) | Initial seeding population. | `> 0`. |
| `t_exp` | `Float64` or `Vector{Float64}` (required) | Expansion phase duration(s). | Scalar or vector; if vector, must align with `Nseed` vector usage. |
| `tmax` | `Float64` or `Vector{Float64}` (required) | Final experiment time horizon(s). | Scalar must be `> 0`; vector values must all be `> 0`, vector length must match effective `n_rep`, and vector mode requires `t_Pass = Float64[]` (no passage events). |
| `t_Pass` | `Float64` or `Vector{Float64}` (required) | Scheduled passage times. | Strictly `> 0`; use `Float64[]` for no passage events. |
| `Nseed` | `Int64` or `Vector{Int64}` (required) | Passage bottleneck sample size(s). | Positive; vector form must align with expansion design. |
| `Nmax` | `Int64` (required) | Population cap trigger for growth/passaging logic. | `> 0` and `Nmax <= Cc`. |
| `Cc` | `Int64` (required) | Carrying-capacity scale in logistic terms. | `> 0`. |
| `treat_ons` | `Vector{Float64}` (required) | Treatment on-times. | Vector of numeric times. |
| `treat_offs` | `Vector{Float64}` (required) | Treatment off-times. | Vector of numeric times. |
| `t_keep` | `Vector{Float64}` (required) | Observation times to retain in compact outputs. | Vector of numeric times. |
| `Nswitch` | `Int64` (required) | Minimum source population for deterministic treatment of each hybrid event. | Positive integer. |
| `Aswitch` | `Float64` (`1000.0`) | Minimum full effective event rate for deterministic treatment, in events per simulation time unit (events/day when rates use days). | Finite and `> 0`. |
| `save_at` | `Float64` (`0.5`) | Hybrid solver save interval. | `> 0`. |
| `n_rep` | `Int64` (`4`) | Number of replicate experiments. | `> 0`; applies per condition in InVivo experiments. |
| `drug_treatment` | `Bool` (`true`) | Whether treatment is applied in main experiment runs. | Boolean. |
| `full_sol` | `Bool` (`false`) | Return full per-time simulated trajectory table. | Boolean. |
| `inc_control` | `Bool` (`false`) | Include untreated controls alongside treated InVivo replicates. | InVivo only; `drug_treatment=false` produces controls only. |
| `inc_pot` | `Bool` (`false`) | Include the common pool in InVivo outputs. | InVivo only; `POT` population rows and, for ABM, `POT_P0` barcode columns. |

Simulation keyword overrides take precedence over stored settings, e.g. `simulate_experiment(model, exp; inc_control=true)`.

### Hybrid event switching

`ResPop` and `ResDmg` classify each division outcome, death, damage, and repair event separately. An event is deterministic only when its source population is at least `Nswitch` **and** its full effective rate is at least `Aswitch`; otherwise it uses a stochastic jump. Full birth/death rates include current treatment, density, fitness, and (in coupled mode) division-outcome adjustments. Independent phenotype conversions use the values described above without birth or density scaling. Biological rates are independent of numerical representation. `SimpleSimParams` defaults to `Nswitch = 100` and `Aswitch = 1000.0`.

`ResPopInVivo` applies the same criteria to each event within each EG stratum, using its EG-specific source population and full effective rate.

Modes are initialized before solving and refreshed at accepted solver steps and after jumps, treatment updates, passaging, and concentration clamps. Check switching thresholds against fully stochastic runs when extinction or founder-arrival timing matters.


## Key API
- `ResPopParams`: ResPop biological parameters.
- `ResDmgParams`: ResDmg biological parameters.
- `ResPopInVivoParams`: ResPop biology with engraftment phenotype parameters.
- `ResPop`: hybrid ODE+jump model.
- `ResDmg`: hybrid ODE+jump model with damaged compartment.
- `ResPop_ABM`: agent-based model.
- `ResDmg_ABM`: agent-based model with damaged compartment.
- `ResPopInVivo`: hybrid model with fixed engraftment phenotypes.
- `ResPopInVivo_ABM`: agent-based model with fixed engraftment phenotypes.
- `engraftment_selection`: sample engraftment survivors from EG-tagged cells or EG-stratified counts.
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
- `src/models/res_pop_in_vivo.jl`: hybrid InVivo model core simulation logic.
- `src/models/res_pop_in_vivo_abm.jl`: ABM InVivo model core simulation logic.
- `src/simulation/simulate.jl`: experiment orchestration drivers and API dispatch.
- `test/runtests.jl`: integration tests.

## Tests
Run tests with: `julia --project -e "using Pkg; Pkg.test()"`

## InVivo Functionality

### Model features

`ResPopInVivo` extends ResPop with a fixed, heritable engraftment phenotype. Each S/R/E phenotype is divided into EG0 and EG1 populations, giving six compartments. EG affects engraftment probability only; growth, drug response, resistant fitness costs and both phenotypic transition modes follow ResPop.

The following model classes use `ResPopInVivoParams`:

- Hybrid class: `ResPopInVivo(params)`
- Agent-Based class: `ResPopInVivo_ABM(params; abm=ABMParams())`

Both classes support `simulate_experiment(model, exp)` and `simulate_simple(model, sim)`. InVivo uses static barcodes in the ABM class; there is no InVivo EvBC variant.

### Expansion, implantation and passaging

Normal experiments begin with scalar or staged expansion. Initial resistant and EG1 cell counts are rounded from `n0*rho` and `n0*fEG1`, and assigned independently. Bottlenecks between expansion stages sample cells without engraftment selection.

After expansion, replicate batches are sampled without replacement from one common pool. With `inc_control=true`, untreated controls (`CO`) are sampled and run before treated replicates (`DT`); `n_rep` applies to each condition. With `drug_treatment=false`, only controls are produced. The pool must contain enough cells for all replicate batches before engraftment.

Engraftment selection is applied after implantation and every actual passage:

- EG1 cells survive independently with probability `pEG`.
- EG0 cells survive independently with probability `pEG*(1-sEG)`.
- Unsuccessful grafts are not replenished; zero survivors represent extinction.

At each passage, cells are sampled to `Nseed` before engraftment selection. Reaching `Nmax` triggers an available passage immediately; otherwise it stops that replicate. Treatment windows remain on continuous experiment time, including across early passages. Vector `tmax` values apply by replicate within each condition and require `t_Pass=Float64[]`.

### Direct splitting after barcoding (ABM only)

Set `split_after_barcoding=true` to skip expansion and split the barcoded starting pool directly into replicate batches. This mode requires an explicit barcode library (`use_lib_probs=true`, with `Nbc` and `bc_probs`).

Resistance and EG1 status are assigned independently to uniformly selected subsets of the unique barcodes observed in the starting pool. `rho` and `fEG1` therefore specify barcode fractions in this mode, and cells sharing a barcode begin with the same resistance and EG status. Replicate sampling and engraftment then proceed as above.

### InVivo outputs and how to access them

`simulate_experiment` returns a `Dict{String, Any}` with population, condition and engraftment outputs:

- `"engraft_df"`: realised surviving counts at implantation (`passage=1`) and subsequent passages, with columns `cond, rep, passage, N_engraft, nEG0_engraft, nEG1_engraft`.
- `"t"`, `"u"`, `"cond"`, `"rep"`: aligned endpoint time, total population, condition and replicate vectors. Each entry describes one condition/replicate endpoint, rather than the observation times in `t_keep`.
- `"sol_df"`: population trajectories with columns `t, nS_EG0, nS_EG1, nR_EG0, nR_EG1, nE_EG0, nE_EG1, nS, nR, nE, n_EG0, n_EG1, N, cond, rep, passage`.

Hybrid experiments include replicate trajectories when `full_sol=true`. With `inc_pot=true`, `sol_df` also contains expansion-pool rows labelled `cond="POT", rep=0, passage=0`; when `full_sol=false`, only those pool rows are included.

ABM experiments additionally return:

- `"lin_df"`: barcode abundance table, with `bc` and condition/replicate/passage columns such as `CO1_P1` and `DT1_P1`.
- `"sub_lin_df"`: optional subsampled barcode counts when `sub_sample_cells=true`, using up to `K` surviving cells per readout.
- `"pheno_bc_df"`: optional phenotype-by-barcode counts when `full_pheno_bc=true`, with columns such as `DT1_P1_S`, `DT1_P1_R`, `DT1_P1_E`. Counts combine EG0 and EG1.

ABM population trajectories are returned unless `just_lin=true`, which also omits the endpoint vectors. `engraft_df` and the requested barcode tables remain available. With `inc_pot=true`, the pool contributes `POT` population rows and `POT_P0` barcode columns; `full_pheno_bc=true` also adds `POT_P0_S`, `POT_P0_R`, `POT_P0_E`. Pool outputs are available with normal expansion or direct splitting.

Simple simulations return `sol_df` with the same population count columns, without experiment condition/replicate labels (and without passage labels for hybrid runs). The ABM class also returns endpoint barcode counts and the optional phenotype/subsample tables; `just_lin=true` omits `sol_df`. Simple simulations do not apply implantation or passage engraftment.

### Minimal InVivo usage example

Using the experimental design `exp` from `Quick start`:

```julia
params_invivo = ResPopInVivoParams(
    b = 0.893,
    d = 0.200,
    rho = 0.1,
    mu = 0.01,
    sig = 0.0,
    del = 0.0,
    al = 0.0,
    Dc = 1.4,
    k = 0.5,
    psi = 0.0,
    drug_effect = :c,
    fEG1 = 0.2,
    pEG = 0.8,
    sEG = 0.5
)

model_invivo = ResPopInVivo_ABM(params_invivo;
    abm = ABMParams(full_pheno_bc=true))
out = simulate_experiment(model_invivo, exp; inc_control=true, inc_pot=true)

out["engraft_df"]
out["pheno_bc_df"]

# The same biological parameters and design can be used for the hybrid class
out_hybrid = simulate_experiment(ResPopInVivo(params_invivo), exp;
    inc_control=true, inc_pot=true, full_sol=true)
```

## Evolvable Barcode (EvBC) Functionality

### Current scope
The EvBC model variants currently provide **ground-truth lineage tracking** for all simulated cells.


### New EvBC model classes
The following Agent-Based model classes are available and dispatched through the same high-level simulation API:

- `ResPop_ABM_EvBC(params; abm=ABMParams())`
- `ResDmg_ABM_EvBC(params; abm=ABMParams())`

You can call either:
- `simulate_experiment(model_evbc, exp)` for full experiment design workflow (expansion/replicates/passaging)
- `simulate_simple(model_evbc, sim)` for a single-run simplified workflow

### EvBC outputs and how to access them

#### 1) Experiment API
`simulate_experiment(model_evbc, exp)` returns a `Dict{String, Any}` with EvBC-specific lineage outputs in addition to standard ABM outputs.

Common keys:
- `"sol_df"`: time-series simulation table
- `"lin_df"`: barcode count table across replicate/passage snapshots
- `"lineage_df"`: ground-truth lineage table

Additional keys:
- `"t"`, `"u"` are included unless `just_lin=true`
- `"sub_lin_df"` is included when `sub_sample_cells=true`

`lineage_df` schema (EvBC):
- `id::Int64`: cell lineage node id, unique across all replicates within one experiment
- `parent_id::Int64`: direct parent id (`0` indicates founder/root record)
- `birth_time::Float64`: simulation time when node was created
- `parent_pheno::String`: parent phenotype at division event
- `child_pheno::String`: child phenotype after transition logic is applied
- `barcode::Float64`: node barcode value (currently inherited ground-truth value)
- `alive_at_end::Bool`: whether that lineage node is alive at simulation end
- `rep::Int64`: replicate index

Passaging preserves retained cells' ids and ancestry; barcodes are unchanged. Lineage tracking starts after the initial expansion and replicate split, so ancestry during expansion is not recorded. Separate simulation calls have independent id spaces.

`sol_df` columns:
- `ResPop_ABM_EvBC`: `t, N, nS, nR, nE, rep`
- `ResDmg_ABM_EvBC`: `t, N, nS, nDS, nDR, nR, rep`

#### 2) Simple API
`simulate_simple(model_evbc, sim)` returns a simplified output dictionary:
- `"sol_df"`
- `"lin_df"`
- `"lineage_df"`

Notes:
- `simulate_simple` does not return `"t"`/`"u"`
- `lineage_df` still includes `rep` (fixed to `1` in simple mode)

### Working with lineage outputs

Tree, phylogeny and parent-barcode utilities reject duplicate cell ids. For older experiment outputs that restarted ids in each replicate, select one replicate before calling these utilities. New experiment outputs have unique ids across replicates and can be used as a combined table.

The following utility functions can be used directly on `lineage_df`:

- `build_phylogeny(lineage_df; extant_only=false)`
    - Returns edge list: `Vector{Tuple{Int64,Int64}}` as `(parent_id, child_id)`

- `build_tree(lineage_df; extant_only=false)`
    - Returns adjacency dictionary: `Dict{Int64,Vector{Int64}}`

- `lineage_to_newick(lineage_df, root_id; extant_only=false)`
    - Returns Newick string for the specified root

- `population_to_newick(lineage_df, root_id; extant_only=false)`
    - Alias of `lineage_to_newick`

- `lineage_node_metadata(lineage_df)`
    - Returns node-level metadata subset (including optional columns such as `alive_at_end` when present)

- `lineage_edge_barcodes(lineage_df)`
    - Returns parent/child edge table with barcode annotations

`extant_only=true` behavior:
- Keeps all nodes with an extant ancestor at the end of the simulation
- Also keeps any internal ancestors required to connect extant nodes
- Excludes extinct-only terminal branches that are not ancestral to extant nodes


### Minimal EvBC usage example

For a runnable EvBC end-to-end example, see `Quick start` (simple simulation first, then full experiment) and use `ResPop_ABM_EvBC` or `ResDmg_ABM_EvBC` as the model class.

After simulation, access lineage outputs via `out["lineage_df"]`, then call utilities such as `build_phylogeny(...)`, `build_tree(...)`, and `lineage_to_newick(...)` (optionally with `extant_only=true`).
