# Coding Conventions

This document describes the coding conventions currently used in `BarcodeModels.jl`. It is intended for agents adding new models or extending existing functionality.

Read [framework_overview.md](framework_overview.md) for architectural intent and [key_interfaces.md](key_interfaces.md) for current scientific and numerical contracts. These conventions explain how they are reflected in code. Apply the latest user-agreed scope and corrections when an earlier plan differs.

## Layer Boundaries

- Keep biological model logic separate from experiment orchestration.
- Keep experimental design logic separate from biological assumptions.
- Keep simulation methodology reusable across compatible biological models.
- Before adding a function, decide whether it belongs to the biological model layer, model-class layer, experimental design layer, or shared helper layer.
- Preserve existing function locations, script boundaries, workflow order, and explicit phenotype equations. Reuse existing helpers and insert focused local functions as needed. Extracting shared engines, relocating functions, or creating helper scripts requires user agreement; scientific readability takes precedence over reducing duplication.

Common placements:

- Biological parameters, rates, transitions, and validation: `src/types/Parameters.jl` and `src/models/`.
- Hybrid state containers and `ComponentArray` conversion: `src/types/State.jl`.
- Model class hierarchy: `src/types/abstract.jl`.
- Public API dispatch: `src/simulation/simulate.jl`.
- Experiment and simple-simulation workflows: `src/simulation/simulate_hybrid.jl`, `src/simulation/simulate_abm.jl`, and `src/simulation/simulate_abm_evbc.jl`.
- Shared ABM output tables and trajectory recording: `src/simulation/abm_outputs.jl`.
- Reusable mechanics: `src/helpers/`.

## Naming

- Biological model names should be short and mechanism-focused, such as `ResPop` or `ResDmg`.
- Concrete implementation names currently combine biological model and model class, such as `ResPop_ABM` and `ResDmg_ABM_EvBC`.
- Parameter structs use `PascalCase` with a `Params` suffix, such as `ResPopParams`.
- State structs use `PascalCase` with a `State` suffix, such as `ResPopState`.
- Internal workflow helpers use a leading underscore, such as `_simulate_experiment_abm`.
- Mutating functions use a trailing `!`, such as `order_cells!`, `birth_mutate_event!`, and `_record_abm_outputs!`.
- Constants use uppercase names, such as `RESPOP_NS_INDEX` and `RESDMG_DRUG_EFFECTS`.
- Name coupled daughter outcomes explicitly, such as `S_asymmetric_R_birth`; name independent conversions as `S_to_R_transition`. Use `_activity` for full biological activity, `_rate` for executed jump rate, `_ode` for representation flags, and `!` for mutating handlers. Reserve switching terminology in numerical code for representation changes. Existing ABM birth-helper names may remain.

## Constructors And Validation

- Parameter keyword constructors should normalise input types to concrete stored types.
- Convert numeric inputs explicitly with `Float64(...)` or `Int64(...)` where existing structs do so.
- Normalise symbolic options at construction time. Existing examples include `normalize_respop_drug_effect`, `normalize_resdmg_drug_effect`, and `normalize_transition_coupling`.
- Model constructors should validate biological parameters before simulation.
- Hybrid model constructors currently use `validate_model_params`.
- ABM constructors currently use stricter validation with `validate_model_params_strict`.
- Keyword model constructors should forward to the parameter constructor:

```julia
ResPop(; kwargs...) = ResPop(ResPopParams(; kwargs...))
ResPop_ABM(; abm::ABMParams = ABMParams(), kwargs...) = ResPop_ABM(ResPopParams(; kwargs...), abm)
```

## Dispatch Style

- Prefer adding methods for new model types over adding conditionals to existing methods.
- Public API dispatch should remain model-agnostic:
  - `simulate_experiment(model, exp; kwargs...)`
  - `simulate_simple(model, sim; kwargs...)`
- `src/simulation/simulate.jl` should route by model class and concrete model type.
- Model-specific core simulation methods should use existing names:
  - `run_model_core_hybrid(model, state, sim; treat=...)`
  - `run_model_core_abm(model, state, sim; treat=...)`
- Use `_simulate_*` methods for internal workflow implementations called by public dispatch wrappers.

## Parameters And Options

- Keep biological parameters in model-specific structs such as `ResPopParams` and `ResDmgParams`.
- Keep experiment-level parameters in `ExperimentParams` and `SimpleSimParams`.
- Keep model-class-specific simulation settings in structs such as `ABMParams` and `SimParams`.
- Do not add model-specific biological assumptions to `ExperimentParams`.
- Do not add experiment workflow settings to biological parameter structs.
- Preserve current public names such as `drug_treatment`, `treat_ons`, `treat_offs`, `transition_coupling`, `Nswitch`, `Aswitch`, `save_at`, and `full_sol`. The deliberate replacement of `N_trans_switch` has no compatibility aliases.
- Consult [key_interfaces.md](key_interfaces.md) for coupling-dependent parameter units, validation, constructor compatibility, and mode propagation.

## Hybrid Model Conventions

- Hybrid states should have a matching state struct in `src/types/State.jl`.
- Hybrid state structs should have a `to_componentarray` method.
- Keep `ComponentArray` field order aligned with index constants in `src/helpers/common_helpers.jl`.
- Add new index constants when adding new hybrid state layouts.
- Reuse existing shared drug-effect helpers in `src/helpers/ode_helpers.jl` without extracting new cross-model machinery as part of a model extension.
- Use callbacks for treatment schedules, passaging, extinction, population limits, and ODE/jump switching, following the existing `ResPop` and `ResDmg` model files.
- Keep the biological event/rate definitions in model files; keep experiment setup and output assembly in simulation files.
- Preserve the complete-event classification, update, projection, and donor-availability contracts in [key_interfaces.md](key_interfaces.md).

## ABM Conventions

- ABM cell types are mutable structs because simulations update cell state in place.
- ABM state structs wrap cell vectors, for example `ResPopABMState`.
- Maintain the live/dead buffer pattern:
  - dead-cell constructors such as `make_dead_cell`
  - `order_cells!`
  - `extend_with_dead_cells!`
  - `alive_cells`
- Increase or validate `Nbuff` rather than silently writing past available dead slots.
- Reuse generic ABM mechanics already in `src/helpers/abm_helpers.jl`; propose any new extraction separately.
- Keep ResPop/ResDmg standard ABM and EvBC types/setup together in each family's `_abm_types.jl` file, and their events together in `_abm_events.jl`, under `src/models/`.
- Keep model wrapper types, constructors, simulation cores and legacy wrappers in the existing ABM/EvBC model scripts.
- Use `!` for event functions and functions that mutate cells, state, counts, or output vectors.
- Apply the shared buffer and proposal-envelope safeguards in [key_interfaces.md](key_interfaces.md) to both coupling modes and both ABM variants.

## EvBC And Lineage Conventions

- EvBC model variants should preserve standard ABM outputs and add lineage outputs.
- Use `LineageRecord` for lineage records where possible.
- `lineage_df` should preserve the expected public columns:
  - `id`
  - `parent_id`
  - `birth_time`
  - `parent_pheno`
  - `child_pheno`
  - `barcode`
  - `alive_at_end`
- Put generic lineage table/tree operations in `src/helpers/lineage_utils.jl`.
- Root-lineage initialisation is shared in lineage utilities through typed forwarding methods in each family's `_abm_types.jl` file, alongside phenotype labels and cell conversions. Birth-record creation stays in the model-family event file.
- Preserve conversion identity and birth-record semantics described in [key_interfaces.md](key_interfaces.md).

## Output Conventions

- Public simulation results are dictionaries keyed by strings.
- Preserve existing output keys:
  - `sol_df`: time-series population output.
  - `lin_df`: barcode abundance output for ABM classes.
  - `sub_lin_df`: optional subsampled barcode output.
  - `lineage_df`: EvBC lineage output.
  - `t` and `u`: compact experiment summaries where used by `simulate_experiment`.
- `simulate_simple` currently returns compact dictionaries without `t` and `u`.
- DataFrame column names are part of the user-facing API; update tests when changing them and identify needed README changes. Edit the README only when requested or included in the approved scope.
- For invalid hybrid inference cases, preserve the existing sentinel-output behaviour unless deliberately redesigning that interface.

## Error Handling And Assertions

- Use constructor validation for user-facing parameter errors.
- Use explicit `error(...)` messages for runtime invariants that should not fail silently.
- Use `@assert` for internal preconditions already established by the surrounding API.
- Prefer informative error messages that name the offending parameter or invariant.
- Avoid silently clamping, dropping, or renaming user inputs unless existing code already does so for that pathway.

## Style

- Follow existing Julia formatting: four-space indentation, spaces around `=`, and spaces after commas.
- Keep keyword arguments explicit and readable for parameter construction.
- Prefer concrete field types in structs, as existing code does.
- Keep comments short and useful. Use comments to explain non-obvious simulation logic, not obvious assignments.
- Prefer focused functions within existing scripts. Reduced duplication alone does not justify new abstractions or extraction.
- Avoid broad refactors while adding a model. Match the current local pattern first.

## Hot-Path Performance

- Keep hot rate/activity functions type-inferable. Retain the existing `Val` specialization for fixed `ComponentArray` field names in local birth/death builders; do not replace it with runtime symbol lookup without checking inference and allocations.
- For a performance change, warm up representative calls before comparing timings and allocations. Separate compilation/GC effects, implementation overhead, and the number of stochastic events selected by the thresholds. Do not infer a speed improvement from an abstraction change alone.

## Tests And Documentation

- Keep `test/runtests.jl` for the user-maintained core package simulation tests; leave it untouched for development checks unless the user explicitly requests a core-suite change. New update/regression checks belong in standalone development files under `test/`, with matching `.gitignore` entries, and must remain outside the Git repository. Do not add development-test includes to `runtests.jl`.
- Add or update tests when changing public constructors, dispatch, output keys, model state layouts, lineage outputs, or experiment orchestration.
- Integration tests should exercise both `simulate_experiment` and `simulate_simple` when a new model is public.
- ABM tests should use small `Nbuff`, `n0`, `Nmax`, and `tmax` values so they run quickly.
- EvBC tests should check lineage columns and at least one parent-child or extant-cell invariant.
- Use the bounded validation workflow in [../barcode-model-implementation/SKILL.md](../barcode-model-implementation/SKILL.md); distinguish focused results from full-suite results and sampling uncertainty from numerical sensitivity.
- Identify README examples or parameter references affected by a change, but update them only when requested or included in the approved scope.
- Keep relevant guidance in `.agents/skills/` and `.github/agents/` aligned with approved changes; preserve older plans as history with clear supersession notes.
