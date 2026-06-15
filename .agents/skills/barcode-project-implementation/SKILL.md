---
name: barcode-project-implementation
description: Implements approved BarcodeModels synthetic-data project specifications by generating scaffold, scripts and parameter tables with consistent naming and path conventions.
---

# Barcode Project Implementation Skill

Use this skill when implementing a project specification into concrete files.

This skill focuses on translating approved project specs into a working scaffold with script placeholders and initial parameter tables.

## Specifications

Before implementing any project, inspect:

- project_specs/

If multiple specifications exist, ask the user which specification should be implemented.

Do not begin file generation until the target specification is explicitly identified.

Treat the approved project specification as the authoritative behavior contract.

## Implementation Workflow

1. Read the target project specification carefully.
2. Validate required sections from barcode-project-design/project_spec_template.md.
3. Identify closest matching template elements in example_project/.
4. Confirm model-family and parameter usage against README.md when interpretation is ambiguous.
5. Produce an implementation plan before writing files.
6. Ask final clarifying questions.
7. Wait for explicit user confirmation to proceed.
8. Implement scaffold, scripts and initial parameter tables.
9. Perform structural and naming verification checks.
10. Produce an implementation record.

## Implementation Rules

- Follow approved specification strictly.
- Reuse existing conventions whenever possible.
- Keep script responsibilities explicit and stage-aware.
- Keep paths consistent with selected path policy.
- Minimize duplication.
- Keep parameter usage in generated scripts consistent with README.md model descriptions.

Uppercase normalization rule:

- Enforce uppercase project-token substitution in all generated project-root path references.
- If project token is PERSIST, generated references must use projects/Barcode_PERSIST.
- Do not retain mixed-case template references unless specification explicitly requires it.

## Verification Requirements

At minimum verify:

- required directory scaffold exists
- required script files exist by stage
- required parameter tables exist and follow naming patterns
- no mixed-case legacy project-root references remain in generated paths
- stage naming is consistent (M0, M1, M2, ...)

## Expected Outputs

When implementing a project specification, provide:

1. Implementation plan.
2. Files created or modified.
3. Parameter tables created.
4. Verification checks and outcomes.
5. Architectural and workflow decisions.

After successful implementation, create a companion implementation record using implementation_record_template.md.

## Related skills

For structure and naming checks, also consult:

- barcode-project-structure
