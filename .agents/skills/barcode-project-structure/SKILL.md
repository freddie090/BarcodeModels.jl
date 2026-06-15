---
name: barcode-project-structure
description: Defines BarcodeModels synthetic-inference project directory structure, naming conventions and artifact contracts using the example_project directory as the baseline reference.
---

# Barcode Project Structure Skill

Use this skill when designing or building a new synthetic-data project scaffold for BarcodeModels workflows.

This skill documents the required directory structure, naming rules, script layout and output artifact contracts.

Baseline reference:

- example_project/

Always use:

- project_structure_reference.md
- naming_conventions.md
- README.md

## Purpose

Use this skill to ensure new projects follow a consistent layout for:

- model modules used by inference pipelines
- simulation and inference scripts
- simulation and fit outputs
- parameter-table organization

If model family definitions or parameter meaning are unclear while designing structure and scripts, consult README.md before finalizing conventions.

## Project-token normalization rule

Example templates may include mixed-case naming.
For new implementations, enforce uppercase project token usage in project-root references:

- projects/Barcode_{UPPERCASE_PROJECT_TOKEN}
- example: projects/Barcode_PERSIST

## Scope

Current baseline is simulation plus multi-stage inference pipeline scaffolding.
SBI experimental-data scripts can remain placeholders unless explicitly requested.

## Related skills

When collecting requirements or implementing files, also use:

- barcode-project-design
- barcode-project-implementation
