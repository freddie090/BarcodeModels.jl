---
name: barcode-model-implementation
description: This skill focuses on translating model specifications into maintainable Julia code that conforms to the existing BarcodeModels architecture.
---

# Barcode Model Implementation

Use this skill when implementing or modifying code within BarcodeModels.jl.

This skill focuses on translating model specifications into maintainable Julia code that conforms to the existing BarcodeModels architecture.

## Model Specifications

Before implementing any model, inspect the contents of:

model_specs/

If multiple specifications exist, ask the user which specification should be implemented.

Do not begin implementation until the target specification has been explicitly identified.

Treat the model specification document as the authoritative description of the required behaviour.

## Implementation Workflow

1. Read the model specification carefully.
2. Identify the closest existing implementation within the codebase.
3. Determine whether the task requires:
   * extending an existing biological model
   * creating a new biological model
   * modifying a simulation method
   * adding new parameters
   * adding new state variables
   * adding new outputs
4. Classify each required change by architectural layer:
   * biological model core
   * simulation method
   * experiment design/orchestration
   * output schema
   * tests/documentation only
5. Produce an implementation plan before writing code.
6. Implement incrementally, reusing existing abstractions wherever possible.
7. Add or update tests.
8. Verify that existing workflows remain functional.

## Integration Audit Checklist

Before considering a model implementation complete, check every relevant integration point:

* parameter type(s), constructors and validation
* state type(s) and component-array conversion
* model class constructor(s)
* package exports and dispatch routing
* Hybrid core implementation
* ABM core implementation
* `simulate_simple`
* `simulate_experiment`
* output schema assembly:
  * `sol_df`
  * `lin_df`
  * `engraft_df` or other model-specific tables
  * `t`, `u` and any metadata vectors
* README or user-facing documentation
* implementation record in `model_specs/`
* focused and full test coverage

## Implementation Rules

* Follow existing framework architecture.
* Reuse existing abstractions whenever possible.
* Avoid introducing unnecessary new patterns.
* Preserve separation between biological models, simulation methods and experimental designs.
* Prefer consistency and maintainability over cleverness.
* Minimise code duplication.
* Keep biological assumptions local to the biological model layer.
* Keep experiment options, replicate/condition orchestration, sampling schedules and output labelling in the simulation/experiment layer.
* For features implemented in more than one model class, define a shared public output contract and keep Hybrid/ABM outputs aligned unless the model specification explicitly requires a difference.

## Output Schema Requirements

When adding or changing outputs:

* Define expected column names and returned dictionary keys before implementation.
* Check simple and experiment workflows separately.
* Include metadata columns/vectors needed to interpret replicate, passage, condition, treatment or sampling status.
* Add internal consistency tests for derived columns rather than only testing that columns exist.
* Preserve backwards-compatible outputs unless the specification explicitly requires a breaking change.

## Testing Requirements

All new functionality should include appropriate tests.

At minimum verify:

* parameter validation
* model initialisation
* expected simulation behaviour
* output generation
* backwards compatibility with existing workflows

For Julia testing:

* Prefer project-local depot/test scripts when the repository provides them.
* Avoid installing dependencies into base Julia unless the user explicitly requests it.
* Run focused tests first when useful, then run the full suite before final handoff.
* Restart Julia after method signature, struct, dispatch or export changes if observed behaviour appears stale.
* Clearly separate related failures from unrelated pre-existing failures before editing unrelated code.

## Expected Outputs

When implementing new functionality, provide:

1. Implementation plan.
2. Modified files.
3. New files created.
4. Tests added or modified.
5. Architectural justification for major design decisions.
6. Documentation and implementation-record updates when public API, output schema, parameters or examples change.

After successful implementation, create a companion implementation record within the model_specs directory.

The implementation record should have the same base name as the specification and append '_implementation'.

Example:

ResPop_MultiMechanism_v1.md

→

ResPop_MultiMechanism_v1_implementation.md

This document serves as a permanent record of what was actually implemented. A template can be found in 'implementation_record_template.md'.

## Related skills

For all implementation tasks, also consult:

- `barcode-model-codebase`
