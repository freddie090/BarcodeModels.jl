# Project Compatibility Guidelines

New project specifications should reuse existing conventions demonstrated in example_project/ wherever possible.

Prefer extending:

- existing directory organization patterns
- existing stage naming conventions
- existing parameter-table naming conventions
- existing artifact placement and fit-stage folder structure

Avoid introducing new naming schemes when existing patterns satisfy requirements.

Default exclusion:

- Do not include harmonic-related scripts, inputs, outputs, or artifact naming in new project specifications unless the user explicitly requests harmonic functionality.

Project designs should explicitly identify:

- reused conventions
- convention deviations and rationale
- any new components required for the requested scope

Uppercase normalization rule:

- Even when source templates use mixed-case project names, new project specs must normalize root naming and path references to projects/Barcode_{UPPERCASE_PROJECT_TOKEN}.
