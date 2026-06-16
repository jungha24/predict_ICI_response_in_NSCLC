# Version 1 Workflow

Version 1 contains the initial NSCLC workflow from QC through subtype discovery and endpoint modeling.

## Main Steps

1. QC and metadata setup
   - Config: `configs/nsclc_v1.json`
   - Notebook: `notebooks/exploratory/01_qc.ipynb`
   - Source: `src/version1/qc.py`

2. Reference atlas integration and annotation
   - R script: `scripts/r/20260309_step1_reference_atlas.R`
   - Supporting metadata summary: `scripts/r/20260309_metadata_summary.R`

3. Subtype discovery
   - Entry point: `src/version1/run_subtype_pipeline.py`
   - Base config: `configs/version1_subtype_base.yaml`
   - Biological SBP config: `configs/version1_subtype_biological_sbp_v1.yaml`

4. Feature importance and clinical association checks
   - `src/version1/permutation_feature_importance.py`
   - `src/version1/cluster_clinical_association.py`

5. Endpoint modeling
   - Entry point: `src/version1/run_endpoint_modeling.py`
   - Base config: `configs/version1_modeling_base.yaml`
   - Biological SBP config: `configs/version1_modeling_biological_sbp_v1.yaml`

## March 19, 2026 Analysis Plan

Version 1 was organized as the initial end-to-end NSCLC ICI workflow. The March 19 plan defined the analysis around the following aims:

- Build a reference atlas for cell-type annotation.
- Discover baseline immune subtypes.
- Define PD-1/PD-L1 ICI responder or durable-benefit endpoints.
- Incorporate clinical covariates into endpoint modeling.
- Explore supervised extensions using early on-treatment change and pre/post-treatment effects.
- Extract biomarker-oriented and drug-targetable features from each subtype.

## Version 1 Feature Blocks

The retained Version 1 immune feature table combined interpretable patient-level immune summaries rather than only using raw cell proportions.

1. Biological SBP and ILR composition features
   - ILR-transformed cell-type composition was considered for the baseline immune composition layer.
   - Biological two-part SBP features were retained for interpretable log-ratio summaries:
     - `LM_SBP_2part_ilr`: lymphoid-myeloid balance.
     - `CD14_CD16_SBP_2part_ilr`: classic/non-classic monocyte balance.

2. Gene-set/module score features
   - T lymphoid axis: `patient_naive_score_mean`, `patient_cytotoxic_score_mean`, and `patient_naive_minus_cytotoxic_score_mean`.
   - Monocyte axis: `patient_cd14_score_mean`, `patient_cd16_score_mean`, and `patient_cd14_minus_cd16_score_mean`.
   - These module scores were calculated from the normalized RNA assay in the discovery cohort, using the relevant immune compartment for each axis.

3. Cell-type latent features
   - Cell-type PCs were considered in the initial plan.
   - The retained biological SBP modeling config uses selected PC features alongside the biological SBP features where specified.

4. Clinical covariates
   - Age, sex, histology, smoking, ECOG, EGFR status, IO line, previous chemotherapy, and previous targeted therapy.
   - PD-L1 TPS was included as an optional categorical covariate in secondary models.

The immune feature source for the biological SBP branch was:

```text
scripts/r/20260320_version1_0_2_step1_immune_feature.R
```

## Endpoint Modeling Strategy

Version 1 used both binary benefit modeling and survival modeling.

- Binary endpoint: DCB versus NCB using the Hsim definition from the clinical metadata, with DCB `n=27` and NCB `n=38` in the March 19 note.
- Survival endpoints:
  - PFS with censoring derived from either progression event or death event.
  - 6-month restricted PFS with administrative censoring at 183 days.
- Binary model: elastic-net logistic regression.
- Survival model: elastic-net Cox regression.
- Preprocessing was fit within each training fold, including one-hot encoding, zero-variance removal, high-correlation pruning at `>0.95`, VIF pruning at `>10`, and elastic-net shrinkage/selection.
- Cross-validation used nested repeated CV with 5 outer splits, 10 outer repeats, and 4 inner splits.
- Stability selection used bootstrap resampling with `n_bootstrap=500` and selection threshold `0.05` in the retained config.
- The retained biological SBP config used:
  - Logistic alpha grid: `1e-3`, `3e-3`, `1e-2`, `3e-2`, `1e-1`, `3e-1`, `1`, `3`.
  - Logistic l1-ratio grid: `0.25`, `0.5`, `0.75`, `0.9`, `0.95`, `0.99`, `1`.
  - Cox penalizer grid: `1e-3`, `3e-3`, `1e-2`, `3e-2`, `1e-1`, `3e-1`, `1`, `3`.
  - Cox l1-ratio grid: `0.25`, `0.5`, `0.75`, `0.9`, `0.95`, `0.99`, `1`.
  - Selection rule: `best`.

## Final Version 1 Composition

In final form, Version 1 consists of:

- A reference-atlas and annotation setup for the pilot NSCLC cohort.
- A baseline subtype-discovery workflow using patient-level immune summaries, including biological SBP/module-score features and selected latent features.
- Clinical association checks for discovered patient clusters.
- Endpoint modeling workflows comparing clinical-only, immune-only, clinical-plus-immune, and PD-L1-augmented models.
- Archived configs, figures, and selected non-sensitive summary outputs needed to understand or rerun the Version 1 analysis.

## Representative Commands

```bash
PYTHONPATH=src conda run -n nsclc-subtype python src/version1/run_subtype_pipeline.py \
  --config configs/version1_subtype_biological_sbp_v1.yaml
```

```bash
conda run -n nsclc-subtype python src/version1/permutation_feature_importance.py \
  --run-dir data/20260309_pilot/results/version1/subtype_discovery/biological_SBP_v1 \
  --n-permutations 50
```

```bash
PYTHONPATH=src conda run -n nsclc-subtype python src/version1/run_endpoint_modeling.py \
  --config configs/version1_modeling_biological_sbp_v1.yaml
```

## Included Files

```text
workflows/version1/
  configs/       # Version 1 config snapshots
  documents/     # Original server-side version notes
  figures/       # Endpoint modeling and metadata-profile figures
  results/       # Selected non-sensitive summary outputs
```

Patient-level feature tables, merged analysis datasets, and raw/intermediate single-cell objects are excluded.
