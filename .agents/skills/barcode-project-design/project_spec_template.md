# Project Specification

## Specification

{ProjectName}_v#.md (example)

## Date Authored

YYYY-MM-DD

## Authored By

barcode_project_designer

## Summary

Brief description of project goals and intended inference-recovery analysis.

## 1. Project Identity

- Project token: {UPPERCASE_PROJECT_TOKEN}
- Project root directory: projects/Barcode_{UPPERCASE_PROJECT_TOKEN}
- Spec version: v#

## 2. Objective And Scope

- Scientific objective
- Pipeline objective
- In-scope model families and stages
- Out-of-scope items for this version

## 3. Directory Plan

Describe required structure under project root:

- abc_models/
- abc_scripts/
- abc_sim_outputs/
- data/
- misc/

Include model-specific subtrees and script locations.

## 4. Parameter Strategy

### 4.1 Fixed Parameters

Table of parameters fixed across the workflow.

### 4.2 Inferred Parameters By Stage

Table by stage (M0, M1, M2, ...), including which parameters vary.

### 4.3 Parameter Table Files

List required files with naming:

- ABC_Sim_{MODEL}_{STAGE}_Param_Table.csv

### 4.4 Parameter Rationale Mapping

- Map each inferred parameter set to its biological rationale.
- Cross-reference model-family and parameter semantics with README.md.
- Note any intentional departures from README-described defaults.

## 5. Script Inventory

For each planned script:

- File path
- Language
- Purpose
- Inputs
- Outputs
- Upstream dependencies
- Downstream consumers

## 6. Execution Workflow

Ordered run sequence across simulation and inference stages.

Include proceed checkpoints where manual confirmation is required.

## 7. Output Artifact Contracts

Define expected output directories and naming patterns, including:

- ground-truth simulation ids
- fit-stage output folders
- stage-explicit artifact names

## 8. Path Policy

- Default path assumption (HPC-style, matching example_project conventions unless overridden)
- Required normalization: use projects/Barcode_{UPPERCASE_PROJECT_TOKEN} in script path references

## 9. Validation And Verification Plan

Define static checks and dry-run checks builder should perform.

## 10. Open Questions And Risks

List unresolved items explicitly.

## 11. Handoff Notes For Builder

Implementation constraints, non-negotiable conventions and final proceed criteria.
