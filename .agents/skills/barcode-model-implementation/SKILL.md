---
name: barcode-model-implementation
description: Implement approved model specifications, package extensions, and scoped maintenance in BarcodeModels.jl while preserving scientific contracts, familiar code structure, and appropriate validation.
---

# Barcode Model Implementation

Use this skill when implementing or modifying code within BarcodeModels.jl.

This skill focuses on translating model specifications into maintainable Julia code that conforms to the existing BarcodeModels architecture.

## Scope And Source Of Intended Behaviour

For a new biological model or substantial scientific extension, inspect `model_specs/` and `package_code_plans/` and identify the target approved specification or plan from the user's request. Ask only when the target or an unresolved scientific choice cannot be established from the session. Clearly scoped maintenance or repairs may proceed from the user's request without creating a model specification or asking them to select an unrelated one.

Apply the latest user-agreed corrections when an earlier plan differs. Preserve earlier plans as history with a clear successor reference; do not restore superseded mechanics. Local/ignored plans may be absent in another checkout, so report missing information rather than assuming it exists.

Read [../barcode-model-codebase/key_interfaces.md](../barcode-model-codebase/key_interfaces.md) for current coupling, complete-event accounting, hybrid boundaries, and ABM/EvBC invariants. Treat those as existing contracts to preserve unless an approved change deliberately replaces them.

## Implementation Workflow

1. Read the target specification, approved plan, or scoped maintenance request.
2. Identify the closest existing implementation within the codebase.
3. Determine whether the task requires:
   * extending an existing biological model
   * creating a new biological model
   * modifying a simulation method
   * adding new parameters
   * adding new state variables
   * adding new outputs
4. Produce a brief implementation plan and select required acceptance checks before writing code.
5. Implement incrementally, reusing existing abstractions wherever possible.
6. Add or update tests.
7. Verify that existing workflows remain functional.

## Implementation Rules

* Follow existing framework architecture.
* Reuse existing abstractions whenever possible.
* Avoid introducing unnecessary new patterns.
* Preserve separation between biological models, simulation methods and experimental designs.
* Prefer consistency and maintainability over cleverness.
* Preserve existing function locations, script boundaries, workflow order, and explicit phenotype equations. Modify or insert focused local functions in place. Reuse existing helpers; extracting shared engines, relocating functions, or creating helper scripts requires user agreement. Scientific readability takes precedence over reducing duplication.
* Keep biological assumptions local to the biological model layer.

## Testing Requirements

Keep new update/regression checks in standalone development files under `test/` and add their paths to `.gitignore`; development tests must remain outside the Git repository. Leave `test/runtests.jl` untouched unless the user explicitly requests a change to the core package simulation suite. Do not include development tests from that entry point. Load their fixture dependencies explicitly when execution is authorised, and report their coverage separately from the core suite.

All new functionality should include appropriate tests. Select checks that verify the affected contracts; routine maintenance does not automatically require new ensembles or a full scientific revalidation.

For new models or substantial extensions, at minimum verify:

* parameter validation
* model initialisation
* expected simulation behaviour
* output generation
* backwards compatibility with existing workflows

For rate/event or numerical representation changes, preserve meaningful complete-event drift accounting, threshold/donor/projection boundaries, conservation and lineage identity, and analytical/numerical checks where relevant. Compare classes only under matching biological formulas, including `R_real`; distinguish sampling uncertainty from numerical sensitivity and scientific adequacy from execution success.

Use the existing isolated Julia workflow described in [../barcode-model-codebase/repository_map.md](../barcode-model-codebase/repository_map.md). Inspect actual entry-point includes and focused-test fixture dependencies; do not assume every test file is executed by `Pkg.test()` or `runtests.jl`. Record seed, environment/settings, executed files, and completion/exit evidence. Keep validation outputs separate from analysis outputs.

Run selected targeted checks followed by any required suite once. Stop when acceptance checks pass; broaden or repeat only for a new change, failure, or unresolved concern. Report a historical full-suite pass separately from later focused checks, rather than presenting both as validation of the latest state. Preserve owner-run processes and compute budgets.

For performance changes, consult [../barcode-model-codebase/coding_conventions.md](../barcode-model-codebase/coding_conventions.md): preserve hot-path inference and existing fixed-field `Val` specialization, warm up representative calls, and compare timings and allocations before claiming improvement. Documentation-only changes need metadata/reference and diff checks, not Julia simulation runs.

## Documentation

Identify README changes needed for public behaviour, but edit the README only when requested or included in the approved scope. Keep relevant agent/skill guidance aligned with approved changes without appending routine implementation details to user documentation.

## Expected Outputs

When implementing new functionality, provide:

1. Implementation plan.
2. Modified files.
3. New files created.
4. Tests added or modified.
5. Architectural justification for major design decisions.

After implementing a model specification or substantial approved package plan, create or update a companion implementation record within `model_specs/`. Use [implementation_record_template.md](implementation_record_template.md). For scoped maintenance, update a relevant existing record when useful; a new model specification/record is not required for every small repair.

The implementation record should have the same base name as the specification and append '_implementation'.

Example:

ResPop_MultiMechanism_v1.md

→

ResPop_MultiMechanism_v1_implementation.md

This document records what was implemented, later agreed corrections, actual validation coverage, and remaining limitations. Do not imply the latest state passed checks that only ran before a follow-up change.

## Related skills

For all implementation tasks, also consult:

- `barcode-model-codebase` at [../barcode-model-codebase/SKILL.md](../barcode-model-codebase/SKILL.md)
