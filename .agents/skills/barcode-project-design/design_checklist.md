# Project Design Checklist

## Project Identity

- What is the project token?
- Is the token uppercase?
- Confirm root name pattern: projects/Barcode_{UPPERCASE_PROJECT_TOKEN}
- What version label should this spec use?

## Objective

- What scientific or methodological question is the project testing?
- What inference-recovery criteria define success?

## Model Scope

- Which biological model family is in scope?
- Is scope single family or multi-family?
- Which sub-model stages are required (M0, M1, M2, ...)?

## Parameter Strategy

- Which parameters are fixed across all simulations?
- Which parameters are inferred by each fit stage?
- How do stage parameter tables expand from simpler to more complex models?
- Do selected parameters and stage combinations align with model definitions in README.md?
- Is biological rationale for each inferred-parameter set explicitly stated?

## Simulation Plan

- How many ground-truth simulations are required per stage?
- Which stage(s) can generate ground-truth data?
- How are simulation ids assigned?

## Pipeline Scripts

- Which script languages are required (shell, Python, Julia, R)?
- What is the execution order across stages?
- Which scripts are local, which are HPC-oriented?

## Directory and Path Policy

- Confirm use of conventions demonstrated in example_project/ for structure.
- Confirm path policy (default HPC-style assumptions unless overridden).
- Confirm uppercase path-token substitution rule in generated scripts.

## Outputs and Artifacts

- Which simulation artifacts are mandatory?
- Which inference artifacts are mandatory?
- What naming conventions should artifacts follow?

## Handoff Completeness

- Are unresolved decisions listed explicitly?
- Is builder proceed gate requirement documented?
- Is the spec complete enough to implement without guesswork?
