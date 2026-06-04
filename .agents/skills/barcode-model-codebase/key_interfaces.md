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

If a change breaks one of these contracts, update the relevant public API, tests, README, and agent docs deliberately.

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

New broad simulation methodologies should usually add a new abstract subtype in `src/models/abstract.jl` and corresponding dispatch paths. New biological models should usually reuse `HybridModel` or `ABMModel` rather than adding a new model class.

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

## Parameter Object Interface

Parameter objects should be concrete, validated, and explicit about layer ownership.

Biological parameter structs should:

- Store biological rates, probabilities, fitness costs, treatment-response parameters, and model-specific biological assumptions.
- Convert inputs to concrete stored types.
- Normalise enumerated options such as `drug_effect`.
- Have validation methods that check biological constraints.

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
- `src/helpers/abm_helpers.jl` for reusable ABM cell-buffer, barcode, concentration, count, and tracking utilities.
- `src/helpers/lineage_utils.jl` for generic lineage table, tree, and Newick utilities.

Helpers should not silently encode a new biological hypothesis. If a helper only applies to one biological model, keep it in that model file until it is genuinely shared.

## Where New Functionality Belongs

Use this placement guide:

- New cell states, transitions, death modes, birth modes, repair modes, resistance mechanisms, or treatment-response hypotheses belong in a biological model implementation and parameter struct.
- New ODE/jump mechanics, ABM mechanics, solver behaviour, cell-buffer behaviour, or lineage bookkeeping belong in the model-class implementation or helpers.
- New seeding, expansion, passaging, observation, replicate, sequencing, sampling, assay, or measurement workflow belongs in experiment parameter types and simulation orchestration.
- New output summaries belong in the relevant simulation workflow, unless they are generic post-processing utilities.
- New lineage analysis utilities belong in `src/helpers/lineage_utils.jl`.
- New plotting behaviour belongs in `src/plotting/`.

When in doubt, start with the narrowest model-specific implementation. Promote code to a helper only once it is clearly shared across models or model classes.
