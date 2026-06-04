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