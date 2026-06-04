# BarcodeModels Compatibility Guidelines

New models should reuse existing concepts where possible.

Prefer extending:

- existing parameter structures
- existing state types
- existing event logic
- existing output formats

Avoid creating entirely new abstractions when an existing mechanism can be reused.

Model designs should explicitly identify:

- reused components
- modified components
- entirely new components