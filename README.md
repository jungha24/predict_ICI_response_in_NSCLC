# Predict ICI Response in NSCLC

## Project at a glance

Goal: Predict ICI response in NSCLC using clinical metadata and PBMC scRNA-seq-derived immune features.

Cohort: 73 NSCLC patients treated with immune checkpoint inhibitors.

<img width="1277" height="805" alt="image" src="https://github.com/user-attachments/assets/9af7c798-fe27-4cff-82c1-408c6ea125fd" />

This project explores whether systemic immune profiles can be used to define patient-level molecular subgroups associated with ICI response in NSCLC.
The long-term goal is to identify responder-enriched immune subgroups and further resolve non-responder subgroups to nominate potential biomarkers, therapeutic targets, and genetic factors underlying distinct immune states. Given the limited cohort size, the current analysis focuses on identifying a responder-associated subgroup, shown as subgroup A in the upper panel. 

## Repository Layout
```text
  README.md
  configs/                 # Original runnable configuration files
  envs/                    # Conda environment specifications
  src/                     # Python source code, organized by workflow version
  scripts/
    r/                     # R scripts used for QC, annotation, feature extraction, and plotting
    utilities/             # Small helper scripts and environment setup
  notebooks/
    exploratory/           # Lightweight exploratory/review notebooks
  docs/
    data_policy.md         # What is excluded from Git and why
    source_notes/          # Original server notes copied into the repository
  workflows/               
```

## Workflow overview

Input modalities:
- Clinical metadata, including PD-L1 TPS and treatment-related variables
- Patient-level immune features derived from PBMC scRNA-seq

Feature library:
- Cell-type composition
- Cell-state aware score
- Curated gene set (pathway) module scores
- De novo gene set module scores
- Latent axes / PC summaries

Modeling strategy:
- Clinical baseline model
- Single immune-feature add-on scan
- Correlation and VIF-based redundancy pruning
- Nested cross-validation and outer validation
- Exploratory FAMD/KNN/Louvain/UMAP clustering using top outer-validation features

## Honest result summary

Nested cross-validation identified multiple immune features with apparent signal when added to the clinical baseline model. The retained top-30 feature table summarizes feature-level ranking stability and ROC-AUC deltas across the three outer-training searches (see results/feature_search_base_v2_single_feature_outer/outer_search_validation/repeated_top_features_to30.csv), and these features were used for exploratory FAMD/KNN/Louvain/UMAP subgroup analysis.

This should be interpreted separately from the final supervised outer-validation result. In the retained 3-fold outer validation for the binary response endpoint, the single selected immune-feature model in each fold had lower held-out ROC-AUC and AUPRC than the clinical baseline model (see results/feature_search_base_v2_single_feature_outer/outer_search_validation/outer_fold_metrics.csv). This suggests that although some immune features showed training-fold signal, the supervised selected-feature model did not generalize robustly in this small cohort.

The main purpose of this repository is therefore not to claim a validated positive predictor. It is to document a leakage-aware small-cohort ML workflow for translational biomarker discovery, including feature search, stability review, exploratory subgroup analysis, negative validation evidence, and limitations.

Reviewer starting points:
- [Reviewer guide](docs/reviewer_guide.md)
- [Limitations](docs/limitations.md)

## Data availability and what can be reviewed

Due to patient privacy and data-use restrictions, raw clinical metadata, cell-level scRNA-seq outputs, full patient-level feature matrices are not included in this repository.

Instead, this repository is intended to document the analysis design and modeling workflow. Reviewers can inspect:

- the overall project structure and workflow organization;
- configuration files used to define feature sets and modeling stages;
- R scripts for quality control, annotation, feature extraction, and visualization;
- Python source code for preprocessing, feature selection, multicollinearity filtering, and model fitting;
- workflow documentation describing staged feature search, redundancy pruning, and outer validation;

The goal of this repository is therefore not to provide a fully public reanalysis of the protected cohort, but to demonstrate the structure, rationale, and implementation of a leakage-aware small-cohort translational modeling pipeline.
