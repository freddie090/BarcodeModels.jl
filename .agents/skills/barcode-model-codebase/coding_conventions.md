# Coding Conventions

This document describes the coding conventions currently used in `BarcodeModels.jl`. It is intended for agents adding new models or extending existing functionality.

The framework overview is the architectural source of truth. These conventions explain how that architecture is reflected in code.

## Layer Boundaries

- Keep biological model logic separate from experiment orchestration.
- Keep experimental design logic separate from biological assumptions.
- Keep simulation methodology reusable across compatible biological models.
- Before adding a function, decide whether it belongs to the biological model layer, model-class layer, experimental design layer, or shared helper layer.

Common placements:

- Biological parameters, rates, transitions, and validation: `src/types/Parameters.jl` and `src/models/`.
- Hybrid state containers and `ComponentArray` conversion: `src/types/State.jl`.
- Model class hierarchy: `src/models/abstract.jl`.
- Public API dispatch: `src/simulation/simulate.jl`.
- Experiment and simple-simulation workflows: `src/simulation/simulate_hybrid.jl`, `src/simulation/simulate_abm.jl`, and `src/simulation/simulate_abm_evbc.jl`.
- Reusable mechanics: `src/helpers/`.

## Naming

- Biological model names should be short and mechanism-focused, such as `ResPop` or `ResDmg`.
- Concrete implementation names currently combine biological model and model class, such as `ResPop_ABM` and `ResDmg_ABM_EvBC`.
- Parameter structs use `PascalCase` with a `Params` suffix, such as `ResPopParams`.
- State structs use `PascalCase` with a `State` suffix, such as `ResPopState`.
- Internal workflow helpers use a leading underscore, such as `_simulate_experiment_abm`.
- Mutating functions use a trailing `!`, such as `order_cells!`, `birth_mutate_event!`, and `_record_abm_outputs!`.
- Constants use uppercase names, such as `RESPOP_NS_INDEX` and `RESDMG_DRUG_EFFECTS`.

## Constructors And Validation

- Parameter keyword constructors should normalise input types to concrete stored types.
- Convert numeric inputs explicitly with `Float64(...)` or `Int64(...)` where existing structs do so.
- Normalise symbolic options at construction time. Existing examples include `normalize_respop_drug_effect` and `normalize_resdmg_drug_effect`.
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
- Preserve existing option names where possible, especially public names such as `drug_treatment`, `treat_ons`, `treat_offs`, `Nswitch`, `N_trans_switch`, `save_at`, and `full_sol`.

## Hybrid Model Conventions

- Hybrid states should have a matching state struct in `src/types/State.jl`.
- Hybrid state structs should have a `to_componentarray` method.
- Keep `ComponentArray` field order aligned with index constants in `src/helpers/common_helpers.jl`.
- Add new index constants when adding new hybrid state layouts.
- Keep shared hybrid drug-effect logic in `src/helpers/ode_helpers.jl` when it can be reused.
- Use callbacks for treatment schedules, passaging, extinction, population limits, and ODE/jump switching, following the existing `ResPop` and `ResDmg` model files.
- Keep the biological event/rate definitions in model files; keep experiment setup and output assembly in simulation files.

## ABM Conventions

- ABM cell types are mutable structs because simulations update cell state in place.
- ABM state structs wrap cell vectors, for example `ResPopABMState`.
- Maintain the live/dead buffer pattern:
  - dead-cell constructors such as `make_dead_cell`
  - `order_cells!`
  - `extend_with_dead_cells!`
  - `alive_cells`
- Increase or validate `Nbuff` rather than silently writing past available dead slots.
- Put generic ABM mechanics in `src/helpers/abm_helpers.jl` when shared across models.
- Keep model-specific ABM events in the model implementation file.
- Use `!` for event functions and functions that mutate cells, state, counts, or output vectors.

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
- Keep model-specific lineage creation inside EvBC model files.

## Output Conventions

- Public simulation results are dictionaries keyed by strings.
- Preserve existing output keys:
  - `sol_df`: time-series population output.
  - `lin_df`: barcode abundance output for ABM classes.
  - `sub_lin_df`: optional subsampled barcode output.
  - `lineage_df`: EvBC lineage output.
  - `t` and `u`: compact experiment summaries where used by `simulate_experiment`.
- `simulate_simple` currently returns compact dictionaries without `t` and `u`.
- DataFrame column names are part of the user-facing API; update tests and README when changing them.
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
- Keep helper functions small when possible, but do not introduce new abstractions unless they reduce real duplication or clarify a framework boundary.
- Avoid broad refactors while adding a model. Match the current local pattern first.

## Tests And Documentation

- Add or update tests when changing public constructors, dispatch, output keys, model state layouts, lineage outputs, or experiment orchestration.
- Integration tests should exercise both `simulate_experiment` and `simulate_simple` when a new model is public.
- ABM tests should use small `Nbuff`, `n0`, `Nmax`, and `tmax` values so they run quickly.
- EvBC tests should check lineage columns and at least one parent-child or extant-cell invariant.
- Update README examples or parameter references when public APIs, output tables, or user-facing parameters change.
- Keep agent-facing docs in `agents_temp/` aligned with code changes when they affect extension workflows.
