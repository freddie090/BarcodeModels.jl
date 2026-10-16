# Repository Map

This document provides a high-level map of the `BarcodeModels.jl` codebase for agents working with the package.

The live repository should be treated as the source of truth. If this document disagrees with the code, follow the code.

For intended changes, use the latest user-agreed scope and corrections; earlier plans remain history. Preserve existing function locations and script boundaries. Reuse existing helpers without treating shared functionality as permission to extract an engine or move workflow code.

## Package Purpose

`BarcodeModels.jl` is a framework for modelling cancer evolution experiments involving cellular barcoding, lineage tracking, treatment response, and resistance evolution.

This map should be read alongside the framework overview. That document defines the core architectural separation between biological hypotheses, simulation methodologies, and experimental designs. This repository map explains where those layers currently live in the codebase and where new functionality should usually be added.

## Framework Layers

The main package layers map to the source tree as follows:

- Biological hypotheses mostly live in `src/types/Parameters.jl`, `src/types/State.jl`, and the biological event/rate definitions inside `src/models/`.
- Simulation methodologies use the hierarchy in `src/types/abstract.jl`, model-class-specific implementations in `src/models/`, and workflow dispatch in `src/simulation/`.
- Experimental designs live in the simulation parameter types in `src/types/Parameters.jl` and the orchestration code in `src/simulation/`.
- Shared utilities live in `src/helpers/`, plotting lives in `src/plotting/`, and verification lives in `test/`.

When adding new functionality, first decide whether it belongs to the biological model, model class, or experimental design layer. Then use the file map below to choose the narrowest place to implement it.

## Quick Repository View

```text
BarcodeModels.jl
|
|-- Project.toml
|   Package metadata and dependencies.
|
|-- README.md
|   User-facing overview, examples, and API notes.
|
|-- src/
|   Package implementation.
|
|   |-- BarcodeModels.jl
|   |   Module entry point, includes, imports, and exports.
|   |
|   |-- types/
|   |   Shared model hierarchy, parameter, state, and output record types.
|   |
|   |   |-- abstract.jl
|   |   |   Common model hierarchy.
|   |   |
|   |   |-- Parameters.jl
|   |   |   Biological, simulation, and experiment parameter structs.
|   |   |
|   |   |-- State.jl
|   |       Hybrid state structs and lineage records.
|   |
|   |-- models/
|   |   Biological model implementations by simulation class.
|   |
|   |   |-- res_pop.jl
|   |   |   Hybrid resistance population model.
|   |   |
|   |   |-- res_dmg.jl
|   |   |   Hybrid resistance-damage model.
|   |   |
|   |   |-- res_pop_abm_types.jl
|   |   |   ResPop standard/EvBC types and setup helpers.
|   |   |
|   |   |-- res_pop_abm_events.jl
|   |   |   ResPop standard/EvBC events.
|   |   |
|   |   |-- res_dmg_abm_types.jl
|   |   |   ResDmg standard/EvBC types and setup helpers.
|   |   |
|   |   |-- res_dmg_abm_events.jl
|   |   |   ResDmg standard/EvBC events.
|   |   |
|   |   |-- res_pop_abm.jl
|   |   |   Agent-based resistance population model.
|   |   |
|   |   |-- res_dmg_abm.jl
|   |   |   Agent-based resistance-damage model.
|   |   |
|   |   |-- res_pop_abm_evbc.jl
|   |   |   Lineage-tracking ResPop ABM.
|   |   |
|   |   |-- res_dmg_abm_evbc.jl
|   |       Lineage-tracking ResDmg ABM.
|   |
|   |-- simulation/
|   |   Public simulation dispatch and experiment orchestration.
|   |
|   |   |-- simulate.jl
|   |   |   Public API dispatch.
|   |   |
|   |   |-- simulate_common.jl
|   |   |   Shared simulation workflow helpers.
|   |   |
|   |   |-- abm_outputs.jl
|   |   |   Shared ABM barcode tables and trajectory output recording.
|   |   |
|   |   |-- simulate_hybrid.jl
|   |   |   Hybrid experiment and simple-simulation workflows.
|   |   |
|   |   |-- simulate_abm.jl
|   |   |   Standard ABM experiment and simple-simulation workflows.
|   |   |
|   |   |-- simulate_abm_evbc.jl
|   |   |   Lineage-aware ABM workflows.
|   |   |
|   |   |-- noise.jl
|   |       Measurement noise helpers.
|   |
|   |-- helpers/
|   |   Reusable cross-model utilities.
|   |
|   |   |-- common_helpers.jl
|   |   |   Indexing, population, logistic, and sampling helpers.
|   |   |
|   |   |-- ode_helpers.jl
|   |   |   Shared hybrid drug-effect helpers.
|   |   |
|   |   |-- abm_helpers.jl
|   |   |   Shared ABM cell-buffer, barcode, and concentration helpers.
|   |   |
|   |   |-- lineage_utils.jl
|   |       Shared root-lineage setup, table, and tree utilities.
|   |
|   |-- plotting/
|       Plotting helpers for simulation outputs.
|
|-- test/
|   Package tests and manual simulation checks.
|
|-- .github/
    Repository automation and GitHub configuration.
```

## Top-Level Layout

- `Project.toml`: package metadata, dependencies, Julia compatibility, and test target.
- `Manifest.toml`: resolved dependency versions for the current environment.
- `README.md`: user-facing package overview, quick starts, parameter reference, and current API documentation.
- `LICENSE.txt`: package licence.
- `src/`: package implementation.
- `test/`: test and manual test entry points.
- `.github/`: repository automation and GitHub configuration.
- `.github/agents/`: model-design and model-coder agent definitions.
- `.agents/skills/`: codebase, model-design, and model-implementation skills and references.
- `package_code_plans/`: dated package-change plans, including complete-event hybrid switching and selectable coupling; some are local/ignored.
- `model_specs/`: dated model specifications and companion implementation records; local/ignored content may be absent in another checkout.
- `run_tests_codex_local.ps1`: existing local isolated Julia test runner, when available.
- `julia_testing_venv/`: local isolated Julia depot/cache, not scientific analysis output.
- `projects/`, `project_specs/`, and `example_project/`: analysis/project workflows and records; remain outside a package-change task unless explicitly included.

## Package Entry Point

- `src/BarcodeModels.jl` defines the `BarcodeModels` module.
- It imports dependencies, includes source files in dependency order, and declares the public exports.
- When adding a new public model, parameter type, state type, helper, or API function, check whether it should be added to the export list here.
- Include order matters because files are loaded by textual inclusion inside one module.

## Core Types

- `src/types/abstract.jl`
  - Defines the abstract model hierarchy:
    - `AbstractBarcodeModel`
    - `HybridModel`
    - `ABMModel`
  - This is the model-class root. Concrete model implementations subtype `HybridModel` or `ABMModel` so the public dispatch layer can route them.
  - Add to this file only when introducing a new broad simulation methodology, not when introducing a new biological model.

- `src/types/Parameters.jl`
  - Defines biological parameter types for current models:
    - `ResPopParams`
    - `ResDmgParams`
  - Defines simulation configuration types for current model classes:
    - `SimParams`
    - `ABMParams`
    - `ABMSimParams`
  - Defines experiment design parameters:
    - `ExperimentParams`
    - `SimpleSimParams`
  - Contains parameter normalisation and validation functions.
  - Biological coupling and simulation thresholds have separate ownership; see [key_interfaces.md](key_interfaces.md) for units, validation, compatibility forwarding, and propagation contracts.
  - Biological model parameters belong here when they describe rates, transition probabilities, fitness costs, treatment response, or other hypothesis-level assumptions.
  - Experimental design parameters belong here when they describe seeding, timing, passaging, treatment windows, observation, or measurement workflow.
  - Model-class settings belong here only when they configure a simulation methodology, such as `ABMParams`.
  - These concepts currently share one file, so keep names and comments clear when extending it.

- `src/types/State.jl`
  - Defines state containers for hybrid models:
    - `ResPopState`
    - `ResDmgState`
  - Defines `LineageRecord` for EvBC lineage output.
  - Provides `to_componentarray` conversions used by hybrid ODE/jump solvers.
  - New hybrid model families usually need a state struct and `to_componentarray` method.

## Model Implementations

- `src/models/res_pop.jl`
  - Hybrid implementation of the `ResPop` biological model.
  - Phenotypes: sensitive (`S`), resistant (`R`), escape (`E`).
  - Contains the concrete `ResPop <: HybridModel` type, parameter validation on construction, component-rate construction, ODE definitions, jump events, callbacks, treatment schedule handling, passaging, and the `run_model_core_hybrid` method.
  - Use this as the closest template for new hybrid compartment models.

- `src/models/res_dmg.jl`
  - Hybrid implementation of the `ResDmg` biological model.
  - Phenotypes/states: sensitive (`S`), damaged-sensitive (`DS`), damaged-resistant (`DR`), resistant (`R`).
  - Contains the concrete `ResDmg <: HybridModel` type, damage and repair rates, ODE definitions, jump events, callbacks, treatment schedule handling, passaging, and the `run_model_core_hybrid` method.

- `src/models/res_pop_abm.jl`
  - Standard agent-based implementation of the `ResPop` biological model.
  - Contains the concrete `ResPop_ABM <: ABMModel` type and constructors, simulation core, `run_model_core_abm`, and the legacy wrapper.
  - Types/setup and events are in the model-family support files below; shared `ABMSimParams` lives in `src/types/Parameters.jl`.
  - Use this as the closest template for new ABM models where cells carry phenotype/barcode identities but not accumulating lineage information used for 'evolving barcodes' (enable tree building).

- `src/models/res_dmg_abm.jl`
  - Standard agent-based implementation of the `ResDmg` biological model.
  - Contains the concrete model type and constructors, simulation core, `run_model_core_abm`, and the legacy wrapper; uses the ResDmg support files below.

- `src/models/res_pop_abm_evbc.jl`
  - EvBC lineage-tracking agent-based implementation of the `ResPop` biological model.
  - Contains the concrete EvBC model type and constructors, lineage-aware simulation core, and `run_model_core_abm`; uses the same ResPop support files as standard ABM.
  - Use this when new functionality needs individual ancestry or cell-history outputs - these are used for 'evolving barcodes' (enable single-cell tree building).

- `src/models/res_dmg_abm_evbc.jl`
  - EvBC lineage-tracking agent-based implementation of the `ResDmg` biological model.
  - Contains the concrete EvBC model type and constructors, lineage-aware simulation core, and `run_model_core_abm`; uses the same ResDmg support files as standard ABM.

- `src/models/res_pop_abm_types.jl` and `src/models/res_dmg_abm_types.jl`
  - Each keeps standard ABM and EvBC sections together for one biological model.
  - Defines cell, count, output, and state types; dead-cell factories; seeding; phenotype labels; typed lineage-initialisation forwarding methods; and cell conversion.
- `src/models/res_pop_abm_events.jl` and `src/models/res_dmg_abm_events.jl`
  - Each contains standard birth, shared phenotype/death handlers, and lineage-aware birth for one biological model; ResDmg also includes damage and repair.
  - Load both family types files before events, and events before the existing ABM model files. Existing model-core order is preserved.

Naming note: existing concrete Julia type names such as `ResPop_ABM` encode both biological model and model class. When designing new concepts, keep the distinction clear in the architecture even if concrete implementation names follow the current package convention.

## Simulation API And Dispatch

- `src/simulation/simulate.jl`
  - Defines the public high-level API:
    - `simulate_experiment(model, exp; kwargs...)`
    - `simulate_simple(model, sim; kwargs...)`
  - Routes by model class:
    - `HybridModel` -> hybrid pipeline
    - `ABMModel` -> ABM pipeline
  - Provides concrete dispatch methods for each existing model.
  - New model classes normally need dispatch entries here.
  - This is the main routing layer between experimental design, biological model, and simulation methodology.

- `src/simulation/simulate_hybrid.jl`
  - Implements experiment and simple-simulation workflows for hybrid models.
  - Calls each model's `run_model_core_hybrid` method.
  - Handles expansion, replicate setup, passage design, treatment flags, compact outputs, and full trajectory outputs for hybrid simulations.
  - Extend this when a hybrid model needs different state initialisation or output table assembly, but keep biological rates/events in `src/models/`.

- `src/simulation/simulate_abm.jl`
  - Implements experiment and simple-simulation workflows for standard ABM models.
  - Handles ABM expansion, replicate splitting, passage loops, barcode summaries, optional subsampling, and output tables.
  - Extend this when standard ABM experiment orchestration or output recording changes.

- `src/simulation/simulate_abm_evbc.jl`
  - Implements experiment and simple-simulation workflows for EvBC ABM models.
  - Adds lineage output construction through `lineage_df` while preserving the ABM-style outputs.
  - Extend this when lineage-aware ABM output structure or passage orchestration changes.

- `src/simulation/simulate_common.jl`
  - Shared simulation helpers.
  - Includes keyword default handling, vector `tmax` validation, replicate-specific time horizons, ABM passage-schedule normalisation, parameter-copy helpers, and drug-effect model cloning.
  - Keep shared scheduling and orchestration helpers here; keep workflow loops in their existing simulation files.

- `src/simulation/abm_outputs.jl`
  - Shared standard ABM and EvBC output helpers: `rep_missings`, `join_dfs`, `get_counts`, both trajectory-append functions, and both model-specific output-recording functions.
  - Included before the standard and EvBC simulation workflows. Passage, expansion, and replicate orchestration remain in those workflow files.

- `src/simulation/noise.jl`
  - Measurement noise helper for simulated outputs.
  - Measurement assumptions belong here or in a future measurement layer, not in biological model files.

## Helper Modules

- `src/helpers/common_helpers.jl`
  - Shared constants and utilities for state indexing, total population calculation, logistic growth factors, passaging draws, and nearest-time lookup.
  - Hybrid models rely on the index constants matching their `ComponentArray` state ordering.
  - Reuse existing mathematical and indexing helpers; propose new extraction separately when required.

- `src/helpers/ode_helpers.jl`
  - Shared hybrid model drug-effect helpers.
  - Defines death-only, birth-only, and combined drug-effect calculations used by hybrid models.
  - Reuse existing helpers. Keep complete-event classification and model-local rate builders in their current model scripts.

- `src/helpers/abm_helpers.jl`
  - Shared ABM utilities for treatment concentration schedules, barcode probability generation and initial sampling, and live/dead cell management.
  - ABM biological events live in their model-family event files; barcode count tables and trajectory-output helpers live in `src/simulation/abm_outputs.jl`.
  - Add ABM mechanics here when they are reusable across ABM biological models.

- `src/helpers/lineage_utils.jl`
  - Utilities for EvBC lineage outputs.
  - Provides shared root-lineage initialisation and `_lineage_df`; family types files retain typed forwarding methods, phenotype labels and cell conversions, while family event files retain birth-record creation.
  - Builds phylogeny/tree structures, Newick strings, node metadata, and edge barcode tables from `lineage_df`.
  - Add lineage-output analysis or conversion helpers here, rather than inside a specific model, when they operate on generic lineage tables.

## Plotting

- `src/plotting/simulation_plots.jl`
  - Plotting helpers for simulation outputs.
  - Public export: `plot_simulation_outputs`.

## Tests

- `test/runtests.jl`
  - User-maintained core package simulation tests used by `Pkg.test()` and the local runner. Leave this file untouched unless the user explicitly requests a core-suite change. Do not include development tests here.

- Other update/regression checks under `test/` are local development tests. Add new development-test paths to `.gitignore` and keep them outside the Git repository. Files already tracked by Git are not untracked merely by adding an ignore entry.
- `test/abm_organization.jl` and `test/lineage_ids.jl`
  - Development checks for seeding, root-lineage metadata, birth ancestry, passage/subsampling outputs and experiment-wide lineage ids.
- `test/drug_effect_boundaries.jl`
  - Development checks for additive drug effects at zero/near-zero baseline rates and ODE/jump drift accounting.
- `test/abm_population_cap.jl`
  - Development checks for immediate initial/cap-reaching ABM stops, phenotype/lineage output and cap-triggered passage timing.
- `test/biological_rate_validation.jl`
  - Development checks for finite/nonnegative biological rates, valid zero boundaries, coupling constraints and hybrid inference sentinels.

- `test/hybrid_events.jl` and `test/hybrid_solver_checks.jl`
  - Focused hybrid accounting, classification, fractional-boundary, analytical, and solver checks; these may be local/ignored.
- `test/transition_coupling.jl`
  - Coupling interfaces, conversion accounting, ABM/EvBC invariants, public APIs, analytical comparisons, and treatment-onset regression checks.
- `test/abm_general_fixes.jl`
  - Baseline validation and shared ordering, buffer, density, treatment, and envelope safeguards.

The core entry point does not include these development files. Do not describe its result as coverage of those files. Inspect their fixture dependencies and load them explicitly for relevant validation when execution is authorised.

- `test/manual_tests.jl`, `test/manual_tests_evbc.jl`, `test/manual_test_helpers.jl`
  - Manual or exploratory test scripts and helpers.
  - These are for use/modification by the user only. Do not add automated tests here. 

Use the existing isolated local workflow when available:

```powershell
./run_tests_codex_local.ps1
```

The runner uses the repository project/Manifest and `julia_testing_venv` depot, disables automatic package precompilation, enables solver output with `BARCODEMODELS_SUPPRESS_SOLVER=false`, and runs with `--compiled-modules=existing`. If the runner is absent, use those same settings and `julia --compiled-modules=existing --project=. test/runtests.jl`. Do not install dependencies or alter the environment merely to validate a documentation change.

`Pkg.test()` is the package entry point in a configured environment, but is not the default isolated local workflow here. See [../barcode-model-implementation/SKILL.md](../barcode-model-implementation/SKILL.md) for selecting checks, recording coverage, and stopping once acceptance checks pass.

## Conventions To Preserve

- Biological models must not contain experiment-specific logic.
- Experimental designs must not contain model-specific biological logic.
- Simulation methodologies should be reusable across biological models.
- Keep model constructors validating biological parameters before simulation.
- Keep public simulation entry points model-agnostic where possible; dispatch should select the implementation.
- Prefer adding methods for new model types over adding conditionals to existing model implementations.
- Keep `ExperimentParams` and `SimpleSimParams` shared unless a model truly needs new experiment-level semantics.
- When adding hybrid states, keep state field order, index constants, and `ComponentArray` conversion aligned.
- When adding ABM states, keep live/dead cell buffer handling consistent with existing helpers.
- Keep output tables compatible with existing user-facing keys such as `sol_df`, `lin_df`, `sub_lin_df`, and EvBC `lineage_df`.
- Identify affected README examples or parameter references, but update the README only when requested or included in the approved scope.
