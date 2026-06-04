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
4. Produce an implementation plan before writing code.
5. Implement incrementally, reusing existing abstractions wherever possible.
6. Add or update tests.
7. Verify that existing workflows remain functional.

## Implementation Rules

* Follow existing framework architecture.
* Reuse existing abstractions whenever possible.
* Avoid introducing unnecessary new patterns.
* Preserve separation between biological models, simulation methods and experimental designs.
* Prefer consistency and maintainability over cleverness.
* Minimise code duplication.
* Keep biological assumptions local to the biological model layer.

## Testing Requirements

All new functionality should include appropriate tests.

At minimum verify:

* parameter validation
* model initialisation
* expected simulation behaviour
* output generation
* backwards compatibility with existing workflows

## Expected Outputs

When implementing new functionality, provide:

1. Implementation plan.
2. Modified files.
3. New files created.
4. Tests added or modified.
5. Architectural justification for major design decisions.

After successful implementation, create a companion implementation record within the model_specs directory.

The implementation record should have the same base name as the specification and append '_implementation'.

Example:

ResPop_MultiMechanism_v1.md

→

ResPop_MultiMechanism_v1_implementation.md

This document serves as a permanent record of what was actually implemented. A template can be found in 'implementation_record_template.md'.

## Related skills

For all implementation tasks, also consult:

- `barcode-models-codebase`