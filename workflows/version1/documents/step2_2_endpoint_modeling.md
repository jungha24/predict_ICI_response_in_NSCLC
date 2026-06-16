# Endpoint Modeling

## Recommended Biological SBP Workflow
```bash
# subtype discovery
PYTHONPATH=src conda run -n nsclc-subtype python src/version1/run_subtype_pipeline.py --config configs/version1_subtype_biological_sbp_v1.yaml
# check importance
conda run -n nsclc-subtype python src/version1/permutation_feature_importance.py \
--run-dir data/20260309_pilot/results/version1/subtype_discovery/biological_SBP_v1 \
--n-permutations 50
# modeling
PYTHONPATH=src conda run -n nsclc-subtype python src/version1/run_endpoint_modeling.py --config configs/version1_modeling_biological_sbp_v1.yaml
```

## Legacy Compositional Workflow
```bash
PYTHONPATH=src conda run -n nsclc-subtype python src/version1/run_endpoint_modeling.py --config configs/version1_modeling_base.yaml
```

## What Changed For `biological_SBP_v1`
- Upstream immune features now come from `Rscripts/20260320_version1_0_2_step1_immune_feature.R`.
- The subtype run uses the updated cell-level metadata file `data/20260309_pilot/results/version1/data_export/20260320_cell_meta.parquet`.
- Instead of recomputing old `prop_*` composition features, the subtype run injects precomputed patient-level biological SBP features from `data/20260309_pilot/results/version1/data_export/20260320_patient_level_immune_features.txt`.
- The biological SBP patient-level features used downstream are:
  - `LM_SBP_2part_ilr`
  - `CD14_CD16_SBP_2part_ilr`
  - `patient_naive_score_mean`
  - `patient_cytotoxic_score_mean`
  - `patient_naive_minus_cytotoxic_score_mean`
  - `patient_cd14_score_mean`
  - `patient_cd16_score_mean`
  - `patient_cd14_minus_cd16_score_mean`
- The subtype output is written to `data/20260309_pilot/results/version1/subtype_discovery/biological_SBP_v1`.
- Endpoint modeling for this variant reads `data/20260309_pilot/results/version1/subtype_discovery/biological_SBP_v1/patient_features.csv`.

## Step By Step
1. `configs/version1_subtype_biological_sbp_v1.yaml`
   - Keeps the existing subtype pipeline structure.
   - Switches `paths.cell_meta` to `20260320_cell_meta.parquet`.
   - Sets `proportion_features.enabled: false`, so the old cell-type proportion block is not generated for this run.
   - Adds `external_patient_features`, which reads the consolidated patient-level biological SBP table and merges those features by `patient_id`.

2. `src/version1/run_subtype_pipeline.py`
   - Still loads counts, genes, barcodes, and cell metadata through the existing subtype pipeline.
   - Builds the discovery cohort exactly as before.
   - Builds patient-level metadata from the selected cells.
   - Builds old proportion features only if `proportion_features.enabled: true`.
   - Loads external patient-level immune features only if `external_patient_features.enabled: true`.
   - Merges patient metadata + external biological SBP features + pseudobulk PC/NMF features into one `patient_features.csv`.
   - Writes `external_patient_features.csv` alongside the standard subtype outputs so the injected biological SBP features are easy to inspect.

3. `configs/version1_modeling_biological_sbp_v1.yaml`
   - Keeps the same clinical input and endpoint definitions as the current endpoint modeling setup.
   - Changes the immune feature table path to `data/20260309_pilot/results/version1/subtype_discovery/biological_SBP_v1/patient_features.csv`.
   - Replaces the old `proportion` block with a new `bio_sbp` block that directly selects the eight biological SBP / module-score features from the subtype output.
   - Keeps the `pc` block, so the endpoint models still compare biological SBP features alone versus biological SBP + pseudobulk PCs. The current config restricts this block to `PC1` and `PC2` via `include_pc_indices: [1, 2]`, and the same block can be changed later with `include_pc_indices` or `max_pc_index`.

4. `src/version1/run_endpoint_modeling.py` and `src/version1/endpoint_modeling/pipeline.py`
   - Read the biological SBP endpoint config.
   - Build one matched patient-level analysis table from:
     - clinical metadata: `20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt`
     - immune features: `data/20260309_pilot/results/version1/subtype_discovery/biological_SBP_v1/patient_features.csv`
   - Resolve the requested analysis designs and deduplicate identical designs across `primary` and `secondary` before fitting.
   - Run elastic-net logistic models for `Binarized_response` and elastic-net Cox models for both `PFS` and `PFS_6m_restricted`.
   - Refit the selected model on full data, run bootstrap stability, and write summary outputs.

## Analyses In `configs/version1_modeling_biological_sbp_v1.yaml`
### Primary
- `clinical_only`
- `clinical_plus_bio_sbp`

### Secondary
- `immune_only_bio_sbp_plus_pc`
- `clinical_plus_bio_sbp`
- `clinical_plus_bio_sbp_plus_pc`
- `clinical_plus_bio_sbp_with_pd_l1`
- `clinical_plus_bio_sbp_plus_pc_with_pd_l1`

## Endpoints
- `Binarized_response`
  - Elastic-net logistic regression on the precomputed `Binarized response` clinical column.
- `PFS`
  - Elastic-net Cox regression using `PFS (Days)` and `pfs_event = PD_Event OR Death_Event`.
- `PFS_6m_restricted`
  - Elastic-net Cox regression that keeps continuous time-to-event information only within the first 183 days; follow-up after day 183 is administratively censored at day 183.

## Current Sparsity / Stability Defaults
- `stability.n_bootstrap: 500`
- `stability.selection_threshold: 0.05`
- `stability.selection_scope: all`
- `logistic.alpha_grid: [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]`
- `logistic.l1_ratio_grid: [0.5, 0.75, 0.9, 0.95, 1.0]`
- `cox.penalizer_grid: [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]`
- `cox.l1_ratio_grid: [0.5, 0.75, 0.9, 0.95, 1.0]`
- `logistic.selection_rule: best`
- `cox.selection_rule: best`
- Available selection rules remain `best` and `one_se` (case-insensitive in config).
- `cox.clinical_penalty_factor: 1.0`

The current default now uses `best`, which favors the strongest inner-CV score directly. If you want a more conservative, sparser selection rule again, switch either logistic or Cox to `selection_rule: one_se`.

## Key Outputs
### Subtype (`biological_SBP_v1`)
- `patient_features.csv`
  - patient-level table used downstream, now containing biological SBP features plus any retained PC/NMF features
- `external_patient_features.csv`
  - direct copy of the injected biological SBP features after aligning them to `analysis_id`
- `selected_features.json`
  - now records whether external patient features were enabled and which columns were used

### Endpoint Modeling (`biological_SBP_v1`)
- `merged_analysis_dataset.tsv`
- `feature_dictionary.csv`
- `analysis_registry.csv`
- `endpoint_registry.csv`
- `*__cv_folds.csv`
- `*__stability.csv`
- `model_metrics_summary.csv`
- `full_data_tuning_summary.csv`
- `full_data_coefficients.csv`
- `full_data_ilr_celltype_weights.csv`
- `comparison_summary.csv`
- `run_manifest.json`

## Notes
- The original compositional endpoint workflow remains available in `configs/version1_modeling_base.yaml`.
- The biological SBP workflow does not overwrite the old scripts; it reuses the same runners with a new config path and a small config-driven extension for external patient-level features.
- `full_data_ilr_celltype_weights.csv` is still written by the endpoint pipeline, but it will only be informative for analyses that actually contain ILR-transformed immune blocks.
