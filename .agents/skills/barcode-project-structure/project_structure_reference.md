# Project Structure Reference

This reference defines the expected scaffold for new BarcodeModels synthetic-inference projects.

For model-family definitions and parameter meaning, use README.md as the primary biological reference when interpreting how structure and scripts should map to model/stage choices.

Baseline template source:

- example_project/

## Root directory

The project root directory name must use uppercase project token normalization:

- projects/Barcode_{UPPERCASE_PROJECT_TOKEN}
- example: projects/Barcode_PERSIST

Within that root, mirror the structure patterns demonstrated in the example project.

## Required top-level directories

- abc_models/
  - Julia model modules consumed by pyABC wrappers and simulation runners.
  - Example module naming pattern in template projects:
    - ABC_DMG_Model.jl
    - ABC_POP_Model.jl
  - Convention:
    - One model module per major inference family.
    - Modules expose the model interface expected by pipeline scripts.

- abc_scripts/
  - Inference and orchestration scripts.
  - Required subdirectories:
    - SIMS/
      - Model-specific script bundles.
      - In the example project this contains:
        - RES_DMG/
        - RES_POP/
    - SBI_EXP/
      - Experimental-data inference scripts.
      - Can remain placeholder in simulation-focused v1.
    - MISC/
      - Utility or ad-hoc helper scripts.

  - Typical simulation pipeline scripts per model family in SIMS/{MODEL}/:
    - Run_ABC_{MODEL}_1_sim_array.jl
    - Run_ABC_{MODEL}_2_pop_abc.py
    - Run_ABC_{MODEL}_3a_pop_post_draws.py
    - Run_ABC_{MODEL}_3b_pop_post_draws.jl
    - Run_ABC_{MODEL}_4_lineage_sim.jl
    - Run_ABC_{MODEL}_5_collect_lin_data.R
    - Run_ABC_{MODEL}_6_sbi_posterior.py
    - Run_ABC_{MODEL}_7_plot_posteriors.R
    - Run_ABC_{MODEL}_pop_plot_abc.py

  - Optional replot/diagnostic scripts may also exist in SIMS/{MODEL}/, for example:
    - Re_Plot_ABC_{MODEL}.py
    - Re_Plot_ABC_{MODEL}.sh

  - Shell launch scripts should be placed under:
    - abc_scripts/SIMS/{MODEL}/shell_scripts/
  - Stage launcher naming pattern:
    - Run_ABC_{MODEL}_{STAGE_ID}_{STEP_NAME}.sh

- abc_sim_outputs/
  - Simulated data and fit outputs.
  - Organized by model family and simulation id.

- data/
  - Include as empty directory unless explicitly populated.

- misc/
  - Include as empty directory unless explicitly populated.

## Required root-level parameter tables

Parameter tables should live at project root unless the specification states otherwise.

Naming pattern:

- ABC_Sim_{MODEL}_{STAGE}_Param_Table.csv

Examples:

- ABC_Sim_RES_DMG_M0_Param_Table.csv
- ABC_Sim_RES_DMG_M1_Param_Table.csv
- ABC_Sim_RES_DMG_M2_Param_Table.csv

## abc_scripts expectations

Model-specific simulation and inference scripts should live under:

- abc_scripts/SIMS/{MODEL}/

Script numbering should reflect pipeline order and match shell-launch orchestration.

In multi-language pipelines, keep one primary responsibility per script:

- Julia: simulation arrays, trajectory generation, lineage simulation
- Python: ABC fitting, posterior processing, SBI model fitting
- R: data collation and posterior plotting
- Shell: batch orchestration on HPC or remote systems

SBI_EXP may remain empty in simulation-only implementations, but directory presence should be preserved.

## abc_sim_outputs expectations

Model-family scope:

- abc_sim_outputs/{MODEL}/

Ground-truth simulation ids:

- {GROUNDTRUTH_STAGE}_sim_{N}

Fit-stage subdirectories:

- fit_{FIT_STAGE}
- examples: fit_M0, fit_M1, fit_M2

Expected outputs at simulation root (abc_sim_outputs/{MODEL}/{GROUNDTRUTH_STAGE}_sim_{N}/) include:

- param_df.csv
- hybrid_sol_df.csv
- hybrid_pop_df.csv
- hybrid_pop_traj_plot.jpg
- abm_sol_df.csv
- abm_pop_df.csv
- abm_lin_df.csv
- abm_pop_traj_plot.jpg

Expected outputs within fit-stage folders (fit_{FIT_STAGE}) typically include:

- ABC database:
  - ABC_{MODEL}_{GROUNDTRUTH_STAGE}_{N}_fit_{FIT_STAGE}.db
- ABC posterior tables:
  - ABC_{MODEL}_{GROUNDTRUTH_STAGE}_{N}_fit_{FIT_STAGE}_gen_{G}_pop_posterior.csv
  - ABC_{MODEL}_{GROUNDTRUTH_STAGE}_{N}_fit_{FIT_STAGE}_pop_posterior.csv
- ABC summary plots:
  - ABC_{MODEL}_{GROUNDTRUTH_STAGE}_{N}_fit_{FIT_STAGE}_pop_posterior.jpg
- SBI inputs and outputs:
  - SBI_{MODEL}_sim_{N}_{GROUNDTRUTH_STAGE}_fit_{FIT_STAGE}_harmonic_inputs.csv
  - SBI_{MODEL}_sim_{N}_{GROUNDTRUTH_STAGE}_fit_{FIT_STAGE}_harmonic_inputs.npz
  - SBI_{MODEL}_sim_{N}_{GROUNDTRUTH_STAGE}_fit_{FIT_STAGE}_npe_post.csv
  - SBI_{MODEL}_sim_{N}_{GROUNDTRUTH_STAGE}_fit_{FIT_STAGE}_npe_pairplot.png
- Lineage inference intermediates:
  - lineage_inference_likelihood_round_0.pkl
  - lineage_inference_posterior_round_0.pkl
- Plot support directory:
  - plots/

Generation-specific ABC files may differ by stage and run settings. The specification should define generation counts explicitly.

## Empty directories to preserve

The following top-level directories should be created even if unused initially:

- data/
- misc/

This keeps scaffolds aligned with expected downstream project evolution.

## Path policy

Current default is to align with HPC-style assumptions shown in the example project, unless user explicitly overrides in specification.

Regardless of path policy, generated project-root references must use the uppercase normalized form:

- projects/Barcode_{UPPERCASE_PROJECT_TOKEN}
