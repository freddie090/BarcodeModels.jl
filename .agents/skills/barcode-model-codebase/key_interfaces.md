# Key Interfaces

This document defines the major contracts between layers in `BarcodeModels.jl`.

Read this alongside `framework_overview.md`, `repository_map.md`, `extension_points.md`, and `coding_conventions.md`. The aim here is not to list every function, but to make clear what each layer must provide and what other layers are allowed to assume.

## Layer Contract Summary

Each extension should preserve these contracts:

- A biological model provides parameters, state meaning, event/rate logic, and a concrete model type.
- A simulation method provides the mechanics for running that biological model under an experiment or simple simulation.
- An experimental design provides seeding, treatment, passaging, observation, and measurement settings.
- Output builders provide stable public dictionaries and tables.
- Helpers provide reusable mechanics without owning biological hypotheses or experiment design.

If a change breaks one of these contracts, address the relevant public API, tests, and agent docs deliberately. Identify any README changes needed, but edit the README only when requested or included in the approved scope. Later user-agreed corrections supersede earlier plans; retain the earlier reasoning as history.

## Biological Model Interface

A biological model is the package representation of a biological hypothesis.

It should provide:

- A model-specific parameter struct, usually in `src/types/Parameters.jl`.
- Parameter normalisation and validation functions.
- State definitions when needed, usually in `src/types/State.jl`.
- A concrete model type in `src/models/`.
- Model constructors that validate biological parameters.
- Event, rate, transition, fitness, and treatment-response logic.
- A core simulation method for each supported model class.

Existing examples:

- `ResPopParams`, `ResPopState`, `ResPop`, `ResPop_ABM`, `ResPop_ABM_EvBC`.
- `ResDmgParams`, `ResDmgState`, `ResDmg`, `ResDmg_ABM`, `ResDmg_ABM_EvBC`.

Biological models should not provide:

- Replicate orchestration.
- Experiment passaging workflow.
- User-facing output table assembly.
- Model-class routing.
- Measurement noise workflow, unless the noise is itself part of the biological hypothesis.

## Model Class Interface

A model class describes how a biological model is simulated.

Current model classes are:

- `HybridModel`: deterministic ODE plus stochastic jump-process simulation.
- `ABMModel`: explicit single-cell agent-based simulation.

Concrete model implementations must subtype the relevant model class:

```julia
struct ResPop <: HybridModel
    params::ResPopParams
end

struct ResPop_ABM <: ABMModel
    params::ResPopParams
    abm::ABMParams
end
```

The model class determines public dispatch through:

```julia
simulate_experiment(model, exp; kwargs...)
simulate_simple(model, sim; kwargs...)
```

New broad simulation methodologies should usually add a new abstract subtype in `src/types/abstract.jl` and corresponding dispatch paths. New biological models should usually reuse `HybridModel` or `ABMModel` rather than adding a new model class.

## Hybrid Simulation Interface

A hybrid model implementation should provide:

- A concrete subtype of `HybridModel`.
- A hybrid state struct.
- A `to_componentarray(state)` method.
- Any index constants needed to read/write the state vector.
- A `run_model_core_hybrid(model, state, sim; treat=...)` method.

`run_model_core_hybrid` is expected to:

- Accept a concrete model, concrete state, and `SimParams`.
- Apply treatment according to `sim` and the `treat` keyword.
- Run the core ODE/jump dynamics.
- Return a solver output compatible with the hybrid simulation workflow.

The hybrid experiment workflows in `src/simulation/simulate_hybrid.jl` are responsible for:

- Creating initial state objects.
- Running expansion phases.
- Splitting replicates.
- Applying passaging schedules.
- Assembling `sol_df`, `t`, and `u` outputs.

Do not put full experiment orchestration inside `run_model_core_hybrid`.

### Resistant Fitness-Cost Option

`simulate_simple` and `simulate_experiment` accept `R_real = "l"` for both model families and all Hybrid, ABM, and EvBC classes. Only the strings `"b"`, `"d"`, and `"l"` are valid. `SimParams`, `ABMSimParams`, expansion helpers and legacy simulation wrappers also default to `"l"`; keep the setting through treatment overrides, staged expansion, replicate/passaging runs.

For undamaged R, before treatment and density scaling, `"b"` gives `(b-(b-d)*del, d)`, `"d"` gives `(b, d+(b-d)*del)`, and `"l"` gives `(b*(1-del), d*(1-del))`. Keep these biological equations visible in each model's rate definitions. All give baseline birth-minus-death `(b-d)*(1-del)`; turnover and coupled daughter-production rates differ. E and damaged-state rules retain their existing formulas.

Model constructors check the option-independent `Dc <= b` bound for birth-only treatment. Simulation setup additionally checks `Dc*(1-psi) <= bR` using the selected option (when `psi < 1`), including drug-effect overrides. The resistant bound must not assume `"l"` before the simulation option is supplied. Before this change, ABM/EvBC defaulted to `"b"` and hybrids were fixed to `"l"`.

### Coupling And Complete-Event Accounting

`transition_coupling` is a biological option shared by ResPop and ResDmg across hybrid, standard ABM, and EvBC classes. `Nswitch` and `Aswitch` configure numerical representation independently of that option.

- `:birth_coupled` is the default. `mu`, `sig`, and ResPop's `al` are division-outcome probabilities. An asymmetric birth adds a daughter of another phenotype without removing the mother phenotype. Ordinary births use the remaining outcome probability: `1 - mu`, `1 - sig - al*gamma` for ResPop R, or `1 - sig` for ResDmg R.
- `:independent_rate` converts an existing cell. The same parameters are rates per simulation time unit (per day in current usage); changing mode does not convert their numerical values. Ordinary births use full division activity and produce same-phenotype daughters.

| Independent event | Full biological activity | Net population effect |
|---|---|---|
| S to R, both models | `nS * mu` | S -1, R +1 |
| R to S, both models | `nR * sig` | R -1, S +1 |
| R to E, ResPop | `nR * al * clamp(gamma, 0, 1)` | R -1, E +1 |

Independent conversions have no birth or density multiplier. S/R conversions have no additional treatment dependence; R/E retains the treatment-level factor. E is terminal for these conversions. ResDmg converts only undamaged S/R cells; damaged cells must repair first. Damage, repair, and death remain distinct events with their existing biological formulas, including damaged-cell death after treatment removal.

Each complete event has one source, one full biological activity, and one population effect. Coupled daughter production adds a cell; conversion conserves total population during ordinary event execution. Explicit phenotype ODE equations and jump handlers must implement the same effects. Do not restore the superseded splitting of biological parameters between `*_o` and `*_j` components.

For undamaged cells, the per-cell drug term is `Dc * clamp(gamma, 0, 1)`, attenuated by `1-psi` for resistant/escape phenotypes. It is additive and does not contain a ratio of baseline birth or death rates. In `:d` mode it adds mortality even when baseline death is zero; in `:b` mode it reduces birth under the existing nonnegative-birth constraints; in `:c` mode it reduces birth to zero and adds any remaining suppression to mortality, including for zero-birth cells. Apply the existing density factor after these adjustments. ResDmg damaged-cell mortality remains `d + Dc` with density scaling and no additional exposure term.

### Hybrid Classification And Updates

Classify each event independently using:

```julia
source_population >= Nswitch && full_event_activity >= Aswitch
```

Equality qualifies. Zero activity is inactive. Recipient population does not affect eligibility. Full activity includes the event's current treatment, density, fitness, and division-outcome factors where applicable, and is calculated independently of its representation flag or donor guard. A deterministic flag activates its ODE contribution and disables its jump rate; a stochastic flag does the converse.

`Nswitch` is a positive integer; `Aswitch` is finite and positive, measured in effective events per simulation time unit. `N_trans_switch` was replaced deliberately: do not add compatibility keywords or property aliases. Output spacing is neither the unit of `Aswitch` nor a dedicated classification clock.

Select the coupling mode at setup and construct only its active jump channels (nine ResPop or twelve ResDmg events). Initialise modes before solving; refresh after accepted solver steps and immediately after jumps, treatment changes, passaging, and numerical clamps. Notify the integrator when state or dynamics change and reset applicable jump caches using the existing solver integration. Preserve the local rate builders, explicit equations, callbacks, and solver setup in their existing model scripts.

### Fractional Counts And Donor Availability

- Preserve a compartment's fraction while any individual active deterministic event changes its count. Check positive activity and nonzero population effect, rather than net drift; deterministic effects can cancel. Independent conversion affects both source and recipient for this check.
- If a compartment is below `Nswitch`, noninteger, and has no active deterministic contribution, stochastic-round it: `floor(n) + 1` with probability equal to the fractional part, otherwise `floor(n)`. Integer counts consume no rounding draw.
- Classify all events, determine all projection eligibility before changing any count, project as one batch, then reclassify once if projection occurred. Do not iterate projection to a fixed point or repeatedly round a stable jump-driven compartment.
- Consuming jumps (death, damage, repair, and independent conversion) execute only when their source count is at least one. Below one, their executed rate is zero while deterministic inflow accumulates. Addition-only births retain their full stochastic activity for positive fractional sources. These guards do not change the activity used for classification.

Retain existing passage sampling/integer conversion, total-population-below-one extinction, population limits, and nonnegativity protection. Subcell donor blocking suppresses depletion while inflow accumulates. Projection preserves instantaneous count expectation but can change a trajectory's total population and does not guarantee unbiased endpoints. Source-only classification permits deterministic inflow to rare recipients; actual threshold recrossings can cause further projections. Do not add recipient thresholds, hysteresis, or extra switching controls without an agreed scientific change.

## ABM Simulation Interface

A standard ABM implementation should provide:

- A concrete subtype of `ABMModel`.
- A cell type.
- A state wrapper around the cell vector.
- Count/output structs as needed.
- Seeding logic.
- Model-specific event functions.
- A `run_model_core_abm(model, state, sim; treat=...)` method.

`run_model_core_abm` is expected to:

- Mutate or advance the provided ABM state according to the model dynamics.
- Respect `ABMSimParams`, treatment settings, buffer size assumptions, and stopping conditions.
- Return a model-specific output object containing tracked time, population, phenotype, and passage vectors.

The ABM workflow in `src/simulation/simulate_abm.jl` is responsible for:

- Running expansion.
- Splitting cells into replicates.
- Managing passage loops.
- Recording barcode abundance.
- Assembling `sol_df`, `lin_df`, and optional `sub_lin_df`.

Do not put experiment-level replicate splitting or public output assembly inside core ABM event functions.

Shared `ABMSimParams` is defined in `src/types/Parameters.jl`. Barcode count tables, trajectory appending, and passage-output recording live in `src/simulation/abm_outputs.jl`; ABM workflow loops remain in their existing simulation files.

ResPop and ResDmg each have model-family `_abm_types.jl` and `_abm_events.jl` support files under `src/models/`, with standard ABM and EvBC sections together. Types/setup load before events, which load before the four existing model scripts. Those scripts retain model wrappers, constructors and simulation cores.

### Shared ABM Safeguards

These numerical and buffer contracts apply to both coupling modes in standard ABM and EvBC:

- `order_cells!` preserves every live cell and its metadata in empty, full, and mixed/unsorted buffers; cleanup replaces only dead cells.
- Check and select a dead slot only inside an accepted birth branch. Non-birth events need no spare slot. A birth without space produces an informative buffer error.
- Stop at `Nmax` before drawing a proposal when the initial population already meets/exceeds the cap. A birth that reaches the cap returns immediately after updating cells/counts and recording the post-birth population at that division time. Cap returns retain the exact time in both `tvec` and `fin_t`, allowing passage handling to resume at the correct event time without an extra waiting-time draw.
- Fixed full-run proposal bounds cover all phenotype rates, fitness costs, treatment states, and ResDmg damaged-cell death even without treatment. The same envelope governs exponential waiting times and disjoint event-selection intervals, so at most one biological event executes per proposal.
- Only independent mode adds conversion bounds: `max(mu, sig + al)` for ResPop or `max(mu, sig)` for ResDmg. This sum bounds competing rates; it is not a probability constraint.
- Clamp density-scaled birth/death factors to zero above capacity. Evaluate treatment-dependent activity at the proposed event time. Preserve the existing concentration lookup and onset-zero initialisation, including treatment-off precedence for simultaneous on/off at zero.
- Below `Nmax`, handle a zero proposal envelope by finishing at the requested horizon with unchanged populations. This is defensive core behaviour, not a reason to bypass public biological validation.
- Independent conversion updates phenotype flags and counts only; it preserves total live population, barcode, alive status, and live/dead position vectors. Independent births reuse existing birth helpers with phenotype-outcome probabilities disabled.

## EvBC Lineage Interface

EvBC model variants extend standard ABM contracts with lineage tracking.

An EvBC model should provide:

- A lineage-aware cell type.
- A lineage-aware state object.
- Initial lineage-state construction.
- Lineage-aware birth and phenotype-change events.
- `LineageRecord` creation or equivalent records compatible with `_lineage_df`.
- A `run_model_core_abm` method for the EvBC model type.

EvBC outputs should preserve standard ABM outputs and add:

- `lineage_df`

The public `lineage_df` should include:

- `id`
- `parent_id`
- `birth_time`
- `parent_pheno`
- `child_pheno`
- `barcode`
- `alive_at_end`

Generic operations on lineage tables belong in `src/helpers/lineage_utils.jl`, not inside one biological model.

Shared root initialisation and `_lineage_df` live in lineage utilities. Typed `_init_lineage_state!(cells, t0)` methods, phenotype labels and cell conversions live in each family's `_abm_types.jl` file. The typed methods forward the model-specific phenotype-label function to the shared initialiser. Birth-record creation lives in each family's `_abm_events.jl` file.

`lineage_df` records births and ancestry, not a complete phenotype-change history. Independent conversions preserve cell ID, parent ID, birth time, next-cell ID, and existing lineage records; they create no birth record. Subsequent births record the parent's phenotype at that division.

EvBC cell ids are unique across all replicates within one experiment. Root initialisation accepts a starting id (default 1); each replicate starts at the previous replicate's returned next unused id, including ids consumed by daughters that later died or were removed by passaging. Passage sampling preserves existing cell identity and ancestry. Parent references may cross replicate labels when such ancestry is supplied; replicate membership is not a cell identity. Initial expansion still uses standard ABM cells, so pre-split ancestry is not recorded. Separate simulation calls have independent id spaces. Tree/phylogeny and parent-barcode utilities reject duplicate ids rather than merging ambiguous ancestors; older combined experiment tables must be filtered to one replicate before these operations.

## Parameter Object Interface

Parameter objects should be concrete, validated, and explicit about layer ownership.

Biological parameter structs should:

- Store biological rates, probabilities, fitness costs, treatment-response parameters, and model-specific biological assumptions.
- Convert inputs to concrete stored types.
- Normalise enumerated options such as `drug_effect`.
- Have validation methods that check biological constraints.

For current ResPop/ResDmg models, validate `b > d` for positive baseline growth in all six constructors; treatment-induced contraction remains allowed. Coupled probabilities retain their bounds, including ResPop's `sig + al <= 1`. Independent rates must be finite and nonnegative, may exceed one, and have no probability sum constraint. Preserve stricter ABM validation and existing drug-effect restrictions.

The common biological validators require finite, nonnegative `b`, `d`, `Dc` and `k`; ResDmg also requires finite, nonnegative `ome`, `zet_S` and `zet_R`. Zero rates remain valid subject to `b > d` and the existing drug-mode constraints. `psi` must be finite and at most 1; finite negative values remain available to the hybrid experiment's existing inference-sentinel path. These checks are applied by model construction; biological parameter structs retain their existing normalisation-only constructors.

Accept the two coupling symbols and their exact string equivalents; reject other modes. Preserve the mode through parameter-copy helpers, expansion, replicate setup, drug-effect overrides, and primitive wrappers. Existing old-arity positional biological-parameter constructors forward with birth coupling as the default.

Experiment parameter structs should:

- Store experiment design settings such as seeding, expansion, treatment windows, passaging, observation, replicates, and assay settings.
- Avoid model-specific biological assumptions.
- Be reusable across compatible biological models and model classes.

Model-class parameter structs should:

- Store simulation-method settings such as ABM buffer size, save cadence, barcode library settings, or hybrid switching thresholds.
- Avoid biological assumptions unless they are genuinely part of a specific simulation methodology.

## State Object Interface

State objects define the current simulation state at the boundary between setup code and core simulation code.

Hybrid states should:

- Store phenotype counts, treatment concentration state, and passage number.
- Convert to `ComponentArray` through `to_componentarray`.
- Keep field order aligned with index constants.

ABM states should:

- Wrap a vector of mutable cell objects.
- Preserve the live/dead buffer convention.
- Be passed into `run_model_core_abm`.

Lineage-aware states should:

- Carry cell vectors plus lineage bookkeeping needed for ancestry reconstruction.
- Preserve compatibility with public lineage output builders.

## Output Interface

Public simulation APIs return dictionaries keyed by strings.

Important output keys:

- `sol_df`: population time-series output.
- `lin_df`: barcode abundance output for ABM workflows.
- `sub_lin_df`: optional subsampled barcode abundance output.
- `lineage_df`: EvBC lineage output.
- `t`: compact experiment summary times.
- `u`: compact experiment summary values.

Output contracts:

- Keep key names stable unless intentionally changing the public API.
- Keep DataFrame column names stable when users or tests rely on them.
- `simulate_experiment` may return compact `t` and `u` outputs.
- `simulate_simple` currently returns compact dictionaries without `t` and `u`.
- ABM and EvBC methods should preserve baseline ABM outputs when adding extra outputs.

When adding a new model, tests should verify the expected keys and key DataFrame columns.

## Public Dispatch Interface

The public API is:

```julia
simulate_experiment(model::Union{HybridModel, ABMModel}, exp::ExperimentParams; kwargs...)
simulate_simple(model::Union{HybridModel, ABMModel}, sim::SimpleSimParams; kwargs...)
```

Dispatch contract:

- Public functions route by model class.
- Class-specific wrappers route by concrete model type.
- Concrete wrappers call internal `_simulate_*` implementations.
- Missing concrete methods should fail with an informative "not implemented" error.

For a new concrete model, add the appropriate methods in `src/simulation/simulate.jl`, for example:

```julia
simulate_experiment_hybrid(model::NewModel, exp::ExperimentParams; kwargs...) =
    _simulate_experiment_hybrid(model, exp; kwargs...)

simulate_simple_hybrid(model::NewModel, sim::SimpleSimParams; kwargs...) =
    _simulate_simple_hybrid(model, sim; kwargs...)
```

## Helper Interface

Helpers should be reusable and layer-neutral within their intended scope.

Use:

- `src/helpers/common_helpers.jl` for shared indexing, population, logistic, sampling, and nearest-time utilities.
- `src/helpers/ode_helpers.jl` for reusable hybrid/ODE drug-effect and rate utilities.
- `src/helpers/abm_helpers.jl` for reusable ABM cell-buffer, initial-barcode sampling, and concentration utilities.
- `src/simulation/abm_outputs.jl` for barcode count tables and ABM trajectory/output-recording helpers.
- `src/helpers/lineage_utils.jl` for generic lineage table, tree, and Newick utilities.

Helpers should not silently encode a new biological hypothesis. Reuse existing helpers; keep new model logic local by default. Sharing across models alone does not authorise extracting an engine, relocating functions, or adding helper scripts. Those structural changes require user agreement; preserve familiar workflow order and explicit phenotype equations.

## Where New Functionality Belongs

Use this placement guide:

- New cell states, transitions, death modes, birth modes, repair modes, resistance mechanisms, or treatment-response hypotheses belong in a biological model implementation and parameter struct.
- New ODE/jump mechanics, ABM mechanics, solver behaviour, cell-buffer behaviour, or lineage bookkeeping belong in the model-class implementation or helpers.
- New seeding, expansion, passaging, observation, replicate, sequencing, sampling, assay, or measurement workflow belongs in experiment parameter types and simulation orchestration.
- New output summaries belong in the relevant simulation workflow, unless they are generic post-processing utilities.
- New lineage analysis utilities belong in `src/helpers/lineage_utils.jl`.
- New plotting behaviour belongs in `src/plotting/`.

When in doubt, start with the narrowest model-specific implementation. Propose any extraction separately when a shared responsibility justifies it, and obtain user agreement before changing existing function locations or script boundaries.
