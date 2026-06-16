# Predict ICI Response in NSCLC

This repository documents an analysis workflow for predicting immune checkpoint inhibitor (ICI) response in a pilot NSCLC cohort using clinical metadata and single-cell RNA-seq-derived immune features.

The workflow includes a feature-library-based response modeling pipeline, including staged feature search, outer validation, and exploratory FAMD/KNN/Louvain/UMAP clustering based on top outer-validation features.

<img width="1277" height="805" alt="image" src="https://github.com/user-attachments/assets/9af7c798-fe27-4cff-82c1-408c6ea125fd" />


This analysis only includes the identification of subgroup A, which represents patients responding to anti-PD-(L)1 therapy, as shown in the upper panel.

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
  workflows/              # Stage 2-5 feature search and outer validation
```
