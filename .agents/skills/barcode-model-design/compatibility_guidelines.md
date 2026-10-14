# BarcodeModels Compatibility Guidelines

New models should reuse existing concepts where possible.

Prefer extending:

- existing parameter structures
- existing state types
- existing event logic
- existing output formats

Avoid creating entirely new abstractions when an existing mechanism can be reused.

Preserve existing function locations, script boundaries, workflow order, and explicit phenotype equations. Reuse existing helpers; keep new model logic local by default. Shared functionality alone does not authorise extracting engines, relocating functions, or adding helper scripts. Obtain user agreement for those structural changes.

Use [../barcode-model-codebase/key_interfaces.md](../barcode-model-codebase/key_interfaces.md) for coupling semantics, numerical event representation, parameter propagation, and ABM/EvBC invariants. Explicitly record any agreed departure from these contracts.

Cross-class comparisons require matching biological formulas, including coupling mode, treatment/density effects, and resistant fitness realisation such as `R_real`. Assess observables with sampling uncertainty and numerical approximation limits; do not require identical stochastic trajectories or claim distribution equivalence from execution success.

Model designs should explicitly identify:

- reused components
- modified components
- entirely new components
