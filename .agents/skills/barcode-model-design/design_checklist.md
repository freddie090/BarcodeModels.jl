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