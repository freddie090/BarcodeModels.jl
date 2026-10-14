# Model Design Checklist

## Scientific Question

- What question is the model intended to answer?
- What hypotheses are being compared?

## Cell States

- What cell states exist?
- Are states mutually exclusive?
- Are states reversible?
- Are states heritable?

## Birth and Death

- Does each state have distinct birth rates?
- Does each state have distinct death rates?
- Are rates density dependent?

## Transitions

- What transitions are possible?
- Are transitions directional?
- Are transitions treatment dependent?
- Do phenotype changes occur as division outcomes, independent conversion, or another explicitly defined process?
- For each event, record source, recipient, full effective activity, parameter units, and net population effect. Does it add a daughter or consume an existing source cell?
- Which treatment, density, and fitness factors apply to that event? Distinguish biological activity from representation flags and donor guards.
- What validation applies in each mode? Does changing mode change parameter units without converting values?
- How is mode selection preserved through constructors, parameter copies, expansion, replicates, wrappers, and overrides?

Use [../barcode-model-codebase/key_interfaces.md](../barcode-model-codebase/key_interfaces.md) for current event contracts rather than inferring semantics from older transition names or plans.

## Mutation

- Are mutations tracked?
- Are mutations explicitly represented?
- Infinite sites assumption?

## Treatment

- How does treatment act?
- Birth effect?
- Death effect?
- Transition effect?

## Barcodes

- Static barcodes?
- Evolvable barcodes?
- CRISPR editing?
- Which events preserve or change barcode identity, cell identity, ancestry, and birth metadata?
- Which events create lineage records? Does the requested output require a phenotype history beyond the existing birth/ancestry table?

## Experimental Design

- Number of replicates?
- Number of time points?
- Number of cells sampled?
- Expansion stage pre-replicates?
- Carrying capacity?
- Maximum time, population size, or both?
- Sub-sample cells when measuring population size or barcodes?
- Treatment timings?

## Possible Outputs

- Population sizes
- Clone sizes
- VAFs
- Lineages (Barcode Identities)
- Trees
- Phenotype frequencies

## Model Class Considerations

- ODE/Hybrid
- Agent-based
- Evolvable Barcodes?
- Which biological formulas must match across classes, including `R_real`, treatment, and density settings?
- Which hybrid representation and fractional-count rules are reused? Any departure from the current contracts requires an explicit scientific decision.
- Which numerical approximation limits matter for the intended observables? Separate them from stochastic sampling uncertainty.

## Handoff And Scope

- Identify existing local functions and scripts to extend; preserve explicit phenotype equations and workflow order.
- Does any proposed shared engine, helper script, or function relocation have user agreement?
- Define acceptance checks before implementation, including accounting, conservation/identity, boundaries, public workflows, and analytical/numerical checks where relevant.
- Record latest agreed corrections, excluded project/example work, and README changes to make only when requested or approved in scope.
