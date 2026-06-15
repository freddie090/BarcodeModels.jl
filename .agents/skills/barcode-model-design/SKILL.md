---
name: barcode-model-design
description: This skill helps agents design new models for the BarcodeModels framework.
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
8. Define public output schemas and returned dictionary keys
9. Classify each requirement by architectural layer
10. Produce implementation specification

Always use:
- design_checklist.md
- model_spec_template.md
- compatibility_guidelines.md

Mandatory two-stage protocol

Stage 1: Discovery (Q&A only)
- Use [design_checklist.md]to interview the user.
- Do not create files during this stage.
- Maintain a list of unknowns until all required fields are answered.
- Ensure the user confirms output tables, output columns, returned vectors, replicate labels, condition labels and sampling/bottleneck semantics when relevant.

Stage 2: Specification drafting
- Allowed only after user confirmation: Proceed to spec.
- Then use [model_spec_template.md] exactly.
- Then run compatibility check using compatibility_guidelines.md.
- Then save one dated spec file in model_specs/.

Blocking condition
- If any required discovery field is unanswered, stop and ask questions.
- Never assume biological hypotheses not provided by the user.

The handoff document should follow the design specified in 'model_spec_template.md' exactly and should be complete and not require any additional information to implement the model.

The handoff document should also make implementation boundaries explicit:

- biological model core requirements
- simulation method requirements
- experiment design/orchestration requirements
- output schema requirements
- testing and documentation requirements

For model specifications targeting multiple classes, such as Hybrid and ABM, specify:

- behaviours that must be equivalent across classes
- outputs that must share the same public schema
- acceptable representation-level differences
- derived output consistency checks that should be tested

For treatment arms, controls, passages, sampling, bottlenecks or replicate metadata, specify:

- where the logic belongs architecturally
- how returned results should be ordered
- which output columns or metadata vectors identify condition, replicate, passage or assay origin

model_specs/

Each specification should be stored as a standalone markdown document following the model specification template. It should be dated.

Before creating a new specification, check whether a specification already exists for the proposed model.

## Related skills

When checking compatibility with the existing package, also use:

- `barcode-model-codebase`
