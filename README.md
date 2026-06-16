# Predict ICI Response in NSCLC

This repository documents an analysis workflow for predicting immune checkpoint inhibitor (ICI) response in a pilot NSCLC cohort using clinical metadata and single-cell RNA-seq-derived immune features.

It contains versioned analysis scripts, configuration snapshots, workflow documentation, selected summary outputs, and representative figures. Raw single-cell data, patient-level metadata, and large intermediate matrices are intentionally excluded.

The workflow includes an initial QC/subtype/modeling pipeline (**Version 1**) and a feature-library-based response modeling pipeline (**Version 2**), including staged feature search, outer validation, and exploratory FAMD/KNN/Louvain/UMAP clustering based on top outer-validation features.

<img width="1233" height="766" alt="image" src="https://github.com/user-attachments/assets/ac8f2f61-010d-412e-9996-1a26fc467392" />

## Repository Layout

```text
.
  README.md
  configs/                 # Original runnable configuration files
  envs/                    # Conda environment specifications
  src/                     # Python source code, organized by workflow version
    version1/
    version2/
  scripts/
    r/                     # R scripts used for QC, annotation, feature extraction, and plotting
    utilities/             # Small helper scripts and environment setup
  notebooks/
    exploratory/           # Lightweight exploratory/review notebooks
  docs/
    data_policy.md         # What is excluded from Git and why
    source_notes/          # Original server notes copied into the repository
  workflows/
    version1/              # QC, subtype discovery, biological SBP, endpoint modeling
    version2/              # Stage 2-5 feature search and outer validation
```

## Workflow Versions

### Version 1

Version 1 built the first end-to-end pipeline:

- QC and annotation setup using `configs/nsclc_v1.json`
- Reference atlas integration and cell type annotation through R scripts
- Patient-level subtype discovery from pseudobulk, PCA/NMF, and cell-type proportion features
- Endpoint modeling for binarized response, PFS, and 6-month restricted PFS
- Later biological SBP update using patient-level immune features instead of older composition-only features

Start here:

- `workflows/version1/README.md`
- `src/version1/`
- `workflows/version1/configs/`

### Version 2

Version 2 was the first dedicated feature-library and feature-search workflow.
Its source protocol was `analysis/n73_manual_trial_20260305/version2/protocol_v2.0.1.pages`.

The protocol separates manually curated immune feature construction from the
Python feature-search pipeline:

- Stage 1 feature library generation from patient-level immune features
  including cell composition, potency/dynamics, curated gene-set scores,
  pseudobulk programs/modules, interaction surrogates, and selected latent
  axes such as PCs/FAMD-derived summaries
- Cell-type-specific feature blocks for B, monocyte, NK, CD4 T, CD8 T,
  nonconventional T, and related manually cleaned subsets
- Clinical baseline covariates including age, sex, histology, smoking, ECOG,
  EGFR status, IO line, previous chemo/targeted therapy, and later PD-L1 TPS
- Stage 2 single-feature add-on scans against the clinical baseline
- Stage 3 redundancy pruning with biological family caps, correlation pruning,
  and VIF pruning
- Stage 4 subset search
- Stage 5 best-subset refits, coefficients, and bootstrap stability
- Outer validation to reduce feature-selection bias

Start here:

- `workflows/version2/README.md`
- `src/version2/`
- `workflows/version2/configs/`

## Data Policy

This repository does not include raw single-cell data, large intermediate objects, or patient-level clinical metadata. Excluded examples include:

- `.rds`, `.RData`, `.h5ad`, `.loom`, `.mtx`, `.mtx.gz`
- raw `data/20260309_pilot/` objects
- merged patient-level analysis datasets
- cell-level score matrices and observation tables
- large stdout logs from long-running jobs

See `docs/data_policy.md` for details.

## Environment

The main conda environment files are in `envs/`.

```bash
conda env create -f envs/nsclc-subtype.yml
conda activate nsclc-subtype
```

Some server runs used:

```bash
PYTHONPATH=src conda run -n nsclc-subtype python <entrypoint> --config <config>
```

## Representative Commands

Version 1 biological SBP subtype discovery:

```bash
PYTHONPATH=src conda run -n nsclc-subtype python src/version1/run_subtype_pipeline.py \
  --config configs/version1_subtype_biological_sbp_v1.yaml
```

Version 1 biological SBP endpoint modeling:

```bash
PYTHONPATH=src conda run -n nsclc-subtype python src/version1/run_endpoint_modeling.py \
  --config configs/version1_modeling_biological_sbp_v1.yaml
```

Version 2 feature search:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src \
conda run -n nsclc-subtype python src/version2/run_feature_search_v2.py \
  --config configs/version2_modeling_base_v2_single_feature_outer.yaml
```

## Included Outputs

The repository includes selected small outputs that help review the workflow:

- model registries, manifests, and metric summaries
- selected feature and pruning summaries
- representative figures and diagnostic plots
- original workflow notes copied from the server

These outputs are documentation artifacts, not a complete reproducibility bundle. Full reruns require the private server-side data directory.
