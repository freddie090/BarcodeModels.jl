# Naming Conventions

## Project naming

- Capture a project token that is uppercase.
- Use uppercase token in path references and generated script variables.

Pattern:

- token: {UPPERCASE_PROJECT_TOKEN}
- root directory name: projects/Barcode_{UPPERCASE_PROJECT_TOKEN}

Example:

- token: PERSIST
- root directory: projects/Barcode_PERSIST

## Model and stage naming

- Model families should use uppercase snake tokens in file names when following existing pipeline patterns.
  - examples: RES_DMG, RES_POP

- Stage names should follow M-index convention.
  - examples: M0, M1, M2

## Parameter table naming

Pattern:

- ABC_Sim_{MODEL}_{STAGE}_Param_Table.csv

## Simulation output naming

Ground-truth simulation id pattern:

- {GROUNDTRUTH_STAGE}_sim_{N}

Fit-stage folder pattern:

- fit_{FIT_STAGE}

## Builder rule

When adapting from template files that include mixed-case project names, rewrite all generated path references to use projects/Barcode_{UPPERCASE_PROJECT_TOKEN} unless explicitly instructed otherwise.
