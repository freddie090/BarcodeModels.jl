---
name: barcode-project-design
description: Guides requirement discovery and specification drafting for BarcodeModels synthetic-data project scaffolds and multi-stage inference pipelines.
---

# Barcode Project Design Skill

Use this skill when defining a new project specification for synthetic-data simulation and multi-stage inference workflows.

The goal is not to write code.

The goal is to convert user requirements into a complete, implementation-ready project specification for barcode_project_builder.

Always use:

- design_checklist.md
- project_spec_template.md
- compatibility_guidelines.md
- README.md

Mandatory two-stage protocol

Stage 1: Discovery (Q and A only)

- Use design_checklist.md to interview the user.
- Do not create implementation files during this stage.
- Maintain a list of unknowns until all required fields are answered.

Stage 2: Specification drafting

- Allowed only after user confirmation: Proceed to spec.
- Use project_spec_template.md exactly.
- Run compatibility check using compatibility_guidelines.md.
- Save one dated spec file in project_specs/.

Blocking condition

- If any required discovery field is unanswered, stop and ask questions.
- Never assume model-family choices, stage scope, or parameter inference strategy.
- If model-family/parameter interpretation is unclear, consult README.md before drafting final specification text.

Naming enforcement

- The project token must be uppercase (for example PERSIST).
- Project root path references in the spec must use projects/Barcode_{UPPERCASE_PROJECT_TOKEN}.

Output rules

- Save specifications as project_specs/YYYY-MM-DD_{ProjectName}_v#.md.
- The specification must be complete enough for builder implementation without additional assumptions.

## Related skills

When checking project-layout compatibility and conventions, also use:

- barcode-project-structure
