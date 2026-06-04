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
8. Produce implementation specification

Always use:
- design_checklist.md
- model_spec_template.md
- compatibility_guidelines.md

The handoff document should follow the design specified in 'model_spec_template.md' exactly and should be complete and not require any additional information to implement the model.

model_specs/

Each specification should be stored as a standalone markdown document following the model specification template. It should be dated.

Before creating a new specification, check whether a specification already exists for the proposed model.

## Related skills

When checking compatibility with the existing package, also use:

- `barcode-models-codebase`