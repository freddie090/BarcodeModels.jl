# Repository Map

This document provides a high-level map of the `BarcodeModels.jl` codebase for agents working with the package.

The live repository should be treated as the source of truth. If this document disagrees with the code, follow the code.

## Package Purpose

`BarcodeModels.jl` is a framework for modelling cancer evolution experiments involving cellular barcoding, lineage tracking, treatment response, and resistance evolution.

This map should be read alongside the framework overview. That document defines the core architectural separation between biological hypotheses, simulation methodologies, and experimental designs. This repository map explains where those layers currently live in the codebase and where new functionality should usually be added.

## Framework Layers

The main package layers map to the source tree as follows:

- Biological hypotheses mostly live in `src/types/Parameters.jl`, `src/types/State.jl`, the biological event/rate definitions inside `src/models/`, and model-specific shared biological event helpers in `src/models/shared/`.
- Simulation methodologies live in `src/models/abstract.jl`, model-class-specific implementations in `src/models/`, and workflow dispatch in `src/simulation/`.
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
|   |   Shared parameter, state, and output record types.
|   |
|   |   |-- Parameters.jl
|   |   |   Biological, simulation, and experiment parameter structs.
|   |   |
|   |   |-- State.jl
|   |       Hybrid state structs and lineage records.
|   |
|   |-- models/
|   |   Biological model implementations by simulation class.
|   |
|   |   |-- abstract.jl
|   |   |   Common model hierarchy.
|   |   |
|   |   |-- res_pop.jl
|   |   |   Hybrid resistance population model.
|   |   |
|   |   |-- res_dmg.jl
|   |   |   Hybrid resistance-damage model.
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
|   |   |   Lineage-tracking ResDmg ABM.
|   |   |
|   |   |-- res_pop_in_vivo.jl
|   |   |   Hybrid ResPop in vivo engraftment model.
|   |   |
|   |   |-- res_pop_in_vivo_abm.jl
|   |   |   ABM ResPop in vivo engraftment model.
|   |   |
|   |   |-- shared/
|   |       Model-specific biological event helpers shared by related implementations.
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
|   |   |-- simulate_hybrid.jl
|   |   |   Hybrid experiment and simple-simulation workflows.
|   |   |
|   |   |-- simulate_abm.jl
|   |   |   Standard ABM experiment and simple-simulation workflows.
|   |   |
|   |   |-- simulate_abm_evbc.jl
|   |   |   Lineage-aware ABM workflows.
|   |   |
|   |   |-- abm_outputs.jl
|   |   |   Shared ABM output table assembly helpers.
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
|   |   |   Shared ABM cell, barcode, and tracking helpers.
|   |   |
|   |   |-- lineage_utils.jl
|   |       Generic lineage table and tree utilities.
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

## Package Entry Point

- `src/BarcodeModels.jl` defines the `BarcodeModels` module.
- It imports dependencies, includes source files in dependency order, and declares the public exports.
- When adding a new public model, parameter type, state type, helper, or API function, check whether it should be added to the export list here.
- Include order matters because files are loaded by textual inclusion inside one module.

## Core Types

- `src/models/abstract.jl`
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
  - Defines experiment design parameters:
    - `ExperimentParams`
    - `SimpleSimParams`
  - Contains parameter normalisation and validation functions.
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
  - Defines cell, phenotype count, output, state, and ABM simulation parameter structures used by the ResPop ABM path.
  - Contains the concrete `ResPop_ABM <: ABMModel` type, `CancerCell`, phenotype count/output structs, ABM state, seeding, core stochastic simulation logic, and `run_model_core_abm`.
  - Use this as the closest template for new ABM models where cells carry phenotype/barcode identities but not accumulating lineage information used for 'evolving barcodes' (enable tree building).

- `src/models/res_dmg_abm.jl`
  - Standard agent-based implementation of the `ResDmg` biological model.
  - Defines ResDmg-specific cell/count/output/state structures and core stochastic simulation logic.

- `src/models/res_pop_abm_evbc.jl`
  - EvBC lineage-tracking agent-based implementation of the `ResPop` biological model.
  - Defines lineage-aware ResPop cells and state, conversion from standard ABM cells, lineage-aware birth/mutation events, and `run_model_core_abm`.
  - Use this when new functionality needs individual ancestry or cell-history outputs - these are used for 'evolving barcodes' (enable single-cell tree building).

- `src/models/res_dmg_abm_evbc.jl`
  - EvBC lineage-tracking agent-based implementation of the `ResDmg` biological model.
  - Defines lineage-aware ResDmg cells and state, conversion from standard ResDmg ABM cells, lineage-aware state wrappers, and `run_model_core_abm`.

- `src/models/res_pop_in_vivo.jl`
  - Hybrid implementation of the ResPop in vivo engraftment model.
  - Splits ResPop S/R/E compartments by static EG0/EG1 engraftment labels while reusing ResPop component-rate logic.

- `src/models/res_pop_in_vivo_abm.jl`
  - ABM implementation of the ResPop in vivo engraftment model.
  - Defines `InVivoCancerCell`, in vivo ABM state/output structs, EG-stratified counts, seeding, and core stochastic simulation logic.

- `src/models/shared/res_pop_abm_events.jl`
  - Model-specific ResPop-family ABM birth/mutation and death event helpers shared by related ABM variants.
  - Standard and EvBC ResPop ABMs share the birth/mutation helper; standard, EvBC, and in vivo ResPop ABMs share the death helper.
  - Keep ResPop S/R/E biological event logic here when it is reused across ResPop-family ABM implementations.

- `src/models/shared/res_dmg_abm_events.jl`
  - Model-specific ResDmg-family ABM birth/mutation, lineage-aware birth, damage, repair, and death event helpers shared by standard and EvBC variants.
  - Keep ResDmg biological event logic here when it is reused across ResDmg-family ABM implementations.

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

- `src/simulation/abm_outputs.jl`
  - Shared ABM output helpers for barcode count tables, optional subsampled lineage counts, ResPop/ResDmg trajectory vector accumulation, and in vivo EG-stratified output tables.
  - Put reusable ABM output table assembly here; keep passage-loop control flow in `simulate_abm.jl` or `simulate_abm_evbc.jl`.

- `src/simulation/simulate_common.jl`
  - Shared simulation helpers.
  - Includes keyword default handling, vector `tmax` validation, passage schedule normalization, replicate-specific time horizons, parameter-copy helpers, and drug-effect model cloning.
  - Put cross-pipeline orchestration helpers here only if they are not specific to one biological model or model class.

- `src/simulation/noise.jl`
  - Measurement noise helper for simulated outputs.
  - Measurement assumptions belong here or in a future measurement layer, not in biological model files.

## Helper Modules

- `src/helpers/common_helpers.jl`
  - Shared constants and utilities for state indexing, total population calculation, logistic growth factors, passaging draws, and nearest-time lookup.
  - Hybrid models rely on the index constants matching their `ComponentArray` state ordering.
  - Add small cross-model mathematical or indexing helpers here.

- `src/helpers/ode_helpers.jl`
  - Shared hybrid model drug-effect helpers.
  - Defines death-only, birth-only, and combined drug-effect calculations used by hybrid models.
  - Add hybrid/ODE helper functions here when they are reusable across biological models.

- `src/helpers/abm_helpers.jl`
  - Shared ABM utilities for treatment concentration schedules, barcode probability generation and sampling, barcode-level phenotype assignment, live/dead cell management, count tables, and tracking vectors.
  - Add ABM mechanics here when they are reusable across ABM biological models.
  - Do not put biological phenotype transition events here; shared biological events belong in `src/models/shared/`.

- `src/helpers/lineage_utils.jl`
  - Utilities for EvBC lineage outputs.
  - Builds public lineage DataFrames from `LineageRecord`s, initializes lineage ids for EvBC cells, and builds phylogeny/tree structures, Newick strings, node metadata, and edge barcode tables from `lineage_df`.
  - Add lineage-output analysis or conversion helpers here, rather than inside a specific model, when they operate on generic lineage tables.

## Plotting

- `src/plotting/simulation_plots.jl`
  - Plotting helpers for simulation outputs.
  - Public export: `plot_simulation_outputs`.

## Tests

- `test/runtests.jl`
  - Main package test entry point used by `Pkg.test()`.

- `test/manual_tests.jl`, `test/manual_tests_evbc.jl`, `test/manual_test_helpers.jl`
  - Manual or exploratory test scripts and helpers.
  - These are for use/modification by the user only. Do not add automated tests here. 

Run the package tests with:

```julia
import Pkg
Pkg.test()
```

or from a shell:

```powershell
julia --project -e "using Pkg; Pkg.test()"
```

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
- Update README examples or parameter reference when a public model or user-facing parameter changes.
