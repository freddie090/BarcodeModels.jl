## Adding A New Model: Usual Touch Points

Before adding code, decide which layer the proposed functionality belongs to:

- Biological model: new cell states, rates, transitions, fitness assumptions, treatment-response mechanisms, or resistance hypothesis.
- Model class: new simulation methodology such as ODE-only, jump-only, hybrid, ABM, lineage-aware ABM, or another computational implementation.
- Experimental design: new seeding, treatment, passaging, observation, sequencing, sampling, or measurement-noise behaviour.

Many architectural problems come from putting functionality in the wrong layer. New designs should strengthen the separation between biological model, model class, and experimental design.

For a new hybrid model family:

1. Add parameter and validation logic in `src/types/Parameters.jl`.
2. Add state and `to_componentarray` logic in `src/types/State.jl`.
3. Add a concrete model type in `src/models/`.
4. Implement `run_model_core_hybrid(model, state, sim; treat=...)`.
5. Add include and export entries in `src/BarcodeModels.jl`.
6. Add `simulate_experiment_hybrid` and `simulate_simple_hybrid` dispatch entries in `src/simulation/simulate.jl`.
7. Extend `src/simulation/simulate_hybrid.jl` if the new model needs model-specific state setup or output assembly.
8. Add tests in `test/runtests.jl` or focused test files.

For a new standard ABM model family:

1. Add or reuse parameter logic in `src/types/Parameters.jl`.
2. Define the model, cell type, count/output state, seeding, event functions, and `run_model_core_abm` in `src/models/`.
3. Reuse helpers from `src/helpers/abm_helpers.jl` where possible for layer-neutral ABM mechanics such as barcode sampling, live/dead buffers, and count tables.
4. If biological event functions are reused by multiple implementations of the same biological model family, put them in `src/models/shared/` rather than `src/helpers/abm_helpers.jl`.
5. Add include and export entries in `src/BarcodeModels.jl`.
6. Add `simulate_experiment_abm` and `simulate_simple_abm` dispatch entries in `src/simulation/simulate.jl`.
7. Extend `src/simulation/simulate_abm.jl` if the new model needs custom experiment orchestration; put reusable output assembly in `src/simulation/abm_outputs.jl`.
8. Add tests that cover seeding, event logic, simple simulation, and experiment simulation.

For a new EvBC ABM model family:

1. Start from the matching standard ABM model if one exists.
2. Add lineage-aware cell fields and state.
3. Maintain `LineageRecord` creation during birth and phenotype-change events.
4. Reuse generic lineage DataFrame and root lineage initialization helpers from `src/helpers/lineage_utils.jl`.
5. Extend `src/simulation/simulate_abm_evbc.jl` for lineage-aware orchestration.
6. Check compatibility with `src/helpers/lineage_utils.jl`.
7. Add tests for `lineage_df`, extant-cell flags, parent-child relationships, and public lineage utilities.
