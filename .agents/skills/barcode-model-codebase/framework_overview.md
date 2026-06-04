# BarcodeModels Framework Overview

## Purpose

BarcodeModels.jl is a framework for modelling cancer evolution experiments involving cellular barcoding, lineage tracking, treatment response and resistance evolution.

The framework is designed around strict separation of:

1. Biological hypotheses
2. Simulation methodologies
3. Experimental designs

New additions should preserve this separation wherever possible.

---

# Core Concepts

## Biological Models

A biological model represents a hypothesis about how cells behave and evolve.

Biological models define:

* cell states
* birth processes
* death processes
* mutation processes
* phenotypic switching processes
* treatment responses
* fitness relationships

Biological models do not define how these processes are simulated.

Examples:

* `ResPop`

  * Resistance population model.
  * Represents resistance through predefined phenotypic compartments that exhibit differential survival during treatment.
  * Cells can only acquire resistance through mutation or phenotypic switching.
  * Also possible for resistant cells to exhibit a fitness cost and lose this cost by transitioning to an 'escape' phenotype.
  * This model can therefore either capture simple pre-existing resistance, or more complex models where slower growing phenotypes subsequently give rise to faster growing resistant ('escape') phenotypes.


* `ResDmg`

  * Resistance-damage model.
  * Represents resistance emergence through differential probability of movement into and out of a drug-induced damage state. 
  * Importantly, this means cells can continue to die when treatment is removed.
  * Most relevant for modelling resistance to chemotherapy and DNA-damaging agents.

Future biological models should represent alternative biological hypotheses rather than alternative simulation methods.

---

## Model Classes

A model class defines how a biological model is simulated.

Model classes are computational implementations rather than biological hypotheses.

Examples:

* Hybrid model

  * Combines deterministic (ODE) and stochastic (jump process) simulation approaches.
  * Only keeps track of individual phenotypic compartments.
  * Increases efficiency by switching between deterministic and stochastic simulation based on population size and transition rate thresholds.
  * The cheaper option but unable to provide cell lineage level outputs or outputs that depend on individual cell histories.

* Agent-based model (ABM)

  * Simulates individual cells explicitly.
  * More computationally expensive but can provide lineage-level outputs and outputs that depend on individual cell histories.

The same biological model should be implementable using the multiple model classes.

For example:

```text
ResPop + Hybrid
ResPop + ABM
ResDmg + Hybrid
ResDmg + ABM
```

The biological hypothesis should remain unchanged across model classes. This is important. The two models should lead to the identical behaviour when it comes to the change in phenotypic compartments over time (which both models are able to capture). 

---

## Experimental Design Layer

Experimental design is independent of both biological models and model classes.

Experimental designs define:

* cell seeding numbers
* expansion periods
* treatment schedules
* passaging schedules
* observation times
* sequencing protocols
* barcode sampling procedures
* measurement noise assumptions

Experimental designs should be reusable across biological models and model classes.

For example:

```text
Experiment A
    ├── ResPop Hybrid
    ├── ResPop ABM
    ├── ResDmg Hybrid
    └── ResDmg ABM
```

(Due to the additional behaviour enabled by the agent-based models, some additional experimental parameters must be passed to it as an additional experimental design layer. For example, the agent-based model requires parameters related to barcode assignment which are not relevant for the hybrid model classes.)

The same experiment should be executable using any compatible model.

This modularity is a core design principle of the framework.

---

# Naming Conventions

## Biological Models

Biological models should use concise names describing the biological mechanism.

Examples:

```text
ResPop
ResDmg
ResSwitch
ResPhylo
ResCRISPR
```

Avoid names that encode implementation details.

Poor examples:

```text
ResPopABM
ResPopHybrid
Resistance_Evolution_Model
model12
```

These describe model classes rather than biological models, are too detailed or too vague to be informative.

---

## Model Classes

Model classes should describe simulation methodology.

Examples:

```text
ABM
Hybrid
ODE
Jump
```

Model classes should not encode biological assumptions.

---

## Experimental Designs

Experimental designs should describe the experiment being simulated and are encoded by the parameter object passed to the simulation function.

Examples:

```text
SimpleSimParams
ExperimentParams
```

Experimental designs should not depend on a specific biological model.

---

# Design Rules

## Rule 1

Biological models must not contain experiment-specific logic.

---

## Rule 2

Experimental designs must not contain model-specific logic.

---

## Rule 3

Simulation methodologies must be reusable across biological models.

---

## Rule 4

New biological models should reuse existing framework abstractions whenever possible.

---

## Rule 5

Before creating a new abstraction, determine whether the proposed functionality belongs to:

* biological model
* model class
* experimental design

Many architectural problems arise when functionality is placed in the wrong layer.

---

# Guidance for Agents

When proposing a new model:

1. Identify the biological hypothesis.
2. Determine whether a new biological model is required.
3. Determine whether an existing model class can be reused.
4. Determine whether an existing experimental design can be reused.
5. Preserve separation between biological model, model class and experiment.

New designs should strengthen, rather than weaken, this modular structure.
