---
name: barcode-model-codebase
description: A skill for understanding and contributing to the codebase of the barcode-based cancer evolution models, including its architecture, conventions, and extension points.
---

# Barcode Models Codebase

Use this skill when working with the BarcodeModels.jl codebase.

Examples:

* understanding existing architecture
* designing new biological models
* implementing new models
* reviewing code
* extending simulation methods
* identifying where functionality should live
* checking consistency with framework conventions

Repository:

```text
BarcodeModels.jl
```

Reading Order

When working on architectural or implementation tasks, consult documents in the following order.

* `framework_overview.md`
  * Architectural principles and naming conventions.

* `repository_map.md`
  * Repository structure and file locations.

* `key_interfaces.md`
  * Core interfaces and responsibilities.

* `extension_points.md`
  * How new functionality should be added.

* `coding_conventions.md`
  * Coding style and framework conventions.

---

Guidance for Agents

Before proposing new functionality:

Determine whether the change belongs to a biological model, simulation method or experimental design.
Determine whether an existing abstraction can be reused.
Determine whether the proposed change preserves framework modularity.
Follow existing naming conventions.
Follow existing testing conventions.
Minimise architectural complexity where possible.

When checking implementation compatibility, audit the relevant extension points:

* parameter structs, constructors and validation
* state structs and state-to-component conversion
* model class constructors and package exports
* model core implementation files
* simulation dispatch in `simulate.jl`
* simple simulation workflows
* experiment workflows
* output schema assembly and returned dictionary keys
* lineage/barcode outputs for ABM or EvBC classes
* documentation and implementation records

For functionality shared across Hybrid, ABM or EvBC classes, verify whether the public outputs should be schema-compatible across classes and whether tests should check derived-output consistency, not just the presence of columns.
