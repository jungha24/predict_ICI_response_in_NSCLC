# Predict ICI Response in NSCLC

## Project at a glance

Goal: Predict ICI response in NSCLC using clinical metadata and PBMC scRNA-seq-derived immune features.

Cohort: 73 NSCLC patients treated with immune checkpoint inhibitors.

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

<img width="1277" height="805" alt="image" src="https://github.com/user-attachments/assets/9af7c798-fe27-4cff-82c1-408c6ea125fd" />

This project explores whether systemic immune profiles can be used to define patient-level molecular subgroups associated with ICI response in NSCLC.
The long-term goal is to identify responder-enriched immune subgroups and further resolve non-responder subgroups to nominate potential biomarkers, therapeutic targets, and genetic factors underlying distinct immune states. Given the limited cohort size, the current analysis focuses on identifying a responder-associated subgroup.

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
