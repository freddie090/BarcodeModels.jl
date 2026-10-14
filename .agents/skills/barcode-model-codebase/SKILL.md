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

* [framework_overview.md](framework_overview.md)
  * Architectural principles and naming conventions.

* [repository_map.md](repository_map.md)
  * Repository structure and file locations.

* [key_interfaces.md](key_interfaces.md)
  * Core interfaces, coupling semantics, hybrid event accounting, fractional-population policy, and ABM/EvBC invariants. Consult these contracts when changing rates, events, or numerical representation.

* [extension_points.md](extension_points.md)
  * How new functionality should be added.

* [coding_conventions.md](coding_conventions.md)
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

Preserve existing function locations, script boundaries, workflow order, and explicit phenotype equations. Reuse existing helpers; keep new model logic local by default. Extracting shared engines, relocating functions, or creating helper scripts requires user agreement. Scientific readability takes precedence over reducing duplication.

Use the current user request and latest agreed plan or correction to establish scope. Earlier plans remain history; do not revive superseded mechanisms. Identify README changes that would be useful, but edit the README only when requested or included in the approved scope.
