---
name: barcode-model-design
description: Design new biological models or substantial scientific extensions for BarcodeModels.jl and produce implementation specifications. Routine maintenance uses the implementation skill.
---

# Barcode Model Design Skill

This skill helps agents design new models for the BarcodeModels framework.

The goal is not to write code.

The goal is to convert a biological question into a complete model specification suitable for implementation.

The design process should follow:

1. Clarify scientific objective
2. Define biological assumptions
3. Define state space
4. Define events
5. Define observables
6. Define parameters
7. Assess compatibility with existing BarcodeModels structures
8. Produce implementation specification

Use:
- [design_checklist.md](design_checklist.md)
- [model_spec_template.md](model_spec_template.md)
- [compatibility_guidelines.md](compatibility_guidelines.md)

Consult [../barcode-model-codebase/key_interfaces.md](../barcode-model-codebase/key_interfaces.md) for current coupling, complete-event accounting, hybrid representation, and ABM/EvBC contracts. Explicitly identify any proposed scientific departure from them.

Preserve existing function locations, script boundaries, workflow order, and explicit phenotype equations. Reuse existing helpers, but do not propose extracting shared engines or relocating functions as an incidental part of a model extension. Such structural changes require user agreement; scientific readability takes precedence over reducing duplication.

The handoff document should follow the design specified in 'model_spec_template.md' exactly and should be complete and not require any additional information to implement the model.

model_specs/

Each specification should be stored as a standalone markdown document following the model specification template. It should be dated.

Before creating a new specification, check whether a specification already exists for the proposed model.

Use the latest agreed decisions and corrections when updating a specification. Preserve superseded plans as history with a clear successor reference. Record the intended implementation scope, excluded project/example workflows, and acceptance checks. Identify any README changes needed, but edit the README only when requested or included in the approved scope.

## Related skills

When checking compatibility with the existing package, also use:

- `barcode-model-codebase` at [../barcode-model-codebase/SKILL.md](../barcode-model-codebase/SKILL.md)
