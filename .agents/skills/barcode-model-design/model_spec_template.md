## Specification

ResPop_MultiMechanism_v1.md (example)

## Date Implemented

YYYY-MM-DD

## Implemented By

barcode_model_design (agent name)

## Summary

Brief description of the model design and its intended purpose.


# Model Specification

## Scientific Objective

## Biological Assumptions

## State Variables

## Events

### Birth

### Death

### Switching

### Mutation

### Treatment

### Complete Event Definitions

| Event | Source | Recipient | Full effective activity and modifiers | Parameter units | Net population effect | Consumes source? | Barcode, identity, and lineage effect |
|---|---|---|---|---|---|---|---|

Define each supported mode, including daughter-addition versus existing-cell conversion, treatment/density/fitness factors, and any states excluded from an event. Distinguish full biological activity from executed jump rates and numerical representation. Reference the current codebase contracts and state any agreed departures.

## Parameters

| Parameter | Meaning and units by mode | Default | Validation and compatibility |
|---|---|---|---|

State whether selecting another mode changes units without converting values. Identify propagation through copies, wrappers, expansion, replicate setup, and overrides.

## Outputs

## Computational Strategy

Describe the supported model classes, reused hybrid/buffer contracts, numerical approximation limits, and settings needed for matched cross-class comparisons (including fitness realisation). Preserve existing function locations, script boundaries, workflow order, and explicit phenotype equations. Identify structural changes separately and record user agreement if required.

### Validation and Acceptance

Define the required checks and stopping criteria before implementation. Include relevant event accounting, boundaries, conservation/identity, constructor/propagation compatibility, and public API checks. Where analytical or cross-class comparisons are needed, specify matching biological settings, sampling uncertainty, and separate numerical sensitivity checks.

## Reuse of Existing Components

## New Components Required

## Questions Remaining

## Handoff Notes for Implementation Agent

Identify the authoritative approved plan and later agreed corrections, files/functions to extend, exclusions, and any README work explicitly included in scope. Earlier plans remain history. Do not infer permission to refactor structure from reuse or reduced duplication.
