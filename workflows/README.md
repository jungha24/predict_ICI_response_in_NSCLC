# Version 2 Workflow


## Protocol Summary

The v2 protocol starts with Stage 1 feature-library construction outside the
main Python search pipeline. It combines clinical metadata with patient-level
immune features derived from manually reviewed cell-type subsets and
scRNA-derived summaries.

Feature-library content includes:

- cell-level composition features
- potency/dynamics and CellRank-like transition or priming summaries where
  available
- curated gene-set and pathway scores, including interferon/inflammatory,
  antigen-presentation, B-cell activation/BCR, myeloid migration/phagocytosis,
  NK cytotoxicity, T-cell helper/regulatory/dysfunction, and checkpoint-related
  programs
- pseudobulk program/module features from PCA/NMF or related program discovery
- interaction surrogates and selected latent axes such as PCs or FAMD-derived
  components
- feature quality filters such as minimum patients per feature, non-zero or
  detection thresholds, and removal of low-value or unstable columns

Cell-type blocks called out in the protocol and companion scripts include B
lineage, monocyte, NK, CD4 T, CD8 T, and nonconventional T features. The
patient-level feature tables are merged by sample/patient ID, and aliases are
used as feature-name prefixes to avoid collisions.

Clinical baseline variables include age, sex, histology, smoking, ECOG, EGFR
status, IO line, previous palliative chemotherapy, previous palliative targeted
therapy, and later PD-L1 TPS in the v2.0.1 update.

# version 2 protocol note

1. Stage 1 feature library
    - patient-level feature library를 먼저 만든 뒤 Python search pipeline에 넣는 구조
    - feature type:
        - cell-level composition
        - potency/dynamics 또는 CellRank-like transition/priming surrogate
        - curated gene/pathway scores
        - pseudobulk de novo program/module
        - selected latent axes, e.g. PC summaries
    - 주요 cell-type block:
        - B lineage
        - monocyte
        - NK
        - CD4 T
        - CD8 T
        - nonconventional T
    - QC/filtering:
        - selected patient만 유지
        - min patients per feature, non-zero/detection 기준, unstable column 제거

2. Stage 2: single-feature add-on scan
   - Baseline: clinical model
   - Test: clinical model plus one candidate immune feature
   - Output: delta metric per endpoint
  
3. Stage 3: redundancy pruning
   - Biological family caps
   - Pairwise correlation pruning
   - VIF pruning

4. Stage 4: subset search
   - Exhaustive search for small subsets
   - Beam-search heuristic for larger subsets

5. Stage 5: best-subset refits
   - Full-data coefficients
   - Bootstrap stability
   - Fold-wise coefficient stability

6. Outer validation
   - Repeats feature selection inside outer folds to reduce selection bias

7. Exploratory patient clustering from outer-validation features
   - Summarize the top 30 feature groups from each of the 3 outer-validation
     folds.
   - Use those features to build a patient-level exploratory matrix.
   - Run FAMD on the mixed clinical + immune feature matrix.
   - Build a k-nearest-neighbor graph on selected FAMD dimensions.
   - Run Louvain clustering and visualize the result with UMAP.
   - Overlay cluster, binarized response, RECIST response, 6-month PFS, and
     PD-L1 TPS on the UMAP.

## Endpoints and Baseline

The v2 search focused on binarized response first, with PFS and 6-month
restricted PFS variants represented in earlier v2 outputs. The main baseline is
`clinical_only`; the 2026-04-03 v2.0.1 update introduced PD-L1 TPS as part of
the baseline feature set and adjusted biological family caps.

## Main Code

```text
src/version2/run_feature_search.py
src/version2/run_feature_search_v2.py
src/version2/report_repeated_outer_top_features.py
src/version2/feature_search/
src/version2/feature_search_base_v2/
```

## Key Configs

```text
configs/version2_modeling_smoketest.yaml
configs/version2_modeling_base.yaml
configs/version2_modeling_base_v2.yaml
configs/version2_modeling_base_v2_single_feature_outer.yaml
```

`version2_modeling_base_v2_single_feature_outer.yaml` is the main retained
configuration for the v2.0.1-style search with outer validation.

## Representative Command

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src \
conda run -n nsclc-subtype python src/version2/run_feature_search_v2.py \
  --config configs/version2_modeling_base_v2_single_feature_outer.yaml

conda run -n nsclc-subtype python src/version2/report_repeated_outer_top_features.py \
  --run-dir data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer \
  --top-n 20
```

## Included Files

```text
workflows/
  configs/       # config snapshots
  documents/     # workflow notes
  figures/       # Feature search and exploratory R plots
  results/       # Selected manifests, summaries, top candidates, and pruning traces
```

Important retained result artifacts include:

```text
workflowsresults/feature_search/base_v2_single_feature_outer/
  stage2_baseline_metrics.csv
  stage2_single_feature_scan.csv
  stage3_pruned_candidates.csv
  stage3_redundancy_pairs.csv
  stage3_vif_trace.csv
  stage5_best_subsets.csv
  stage5_best_subset_coefficients.csv
  stage5_best_subset_bootstrap_stability.csv
  stage5_best_subset_coefficient_stability.csv
  outer_search_validation/
    repeated_top_features__top30.csv
    repeated_top_features__top30_details.csv
```

Outer-validation feature and clustering figures include:

```text
workflows/figures/exploratory_rplots/
  20260427_version2_trial1_fold_group_rank.pdf
  20260429_outervalidation_fold_rocauc.pdf
  20260429_version2_trial1_FAMD_FD12.pdf
  20260429_version2_trial1_FAMD_FD23.pdf
  20260429_version2_trial1_FAMD_FD1_explain.pdf
  20260429_version2_trial1_FAMD_FD2_explain.pdf
  20260429_version2_trial1_FAMD_umap.pdf
  20260429_version2_trial1_FAMD_umap_binarized_response.pdf
  20260429_version2_trial1_FAMD_umap_PFS_6mo.pdf
  20260429_version2_trial1_FAMD_umap_RECIST.pdf
  20260429_version2_trial1_FAMD_umap_PDL1_TPS.pdf
  20260429_version2_trial1_FAMD_cluster_DEF*.pdf
```

The clustering script is retained at
`scripts/r/20260408_feature_based_clustering.R`. Its top-30 branch reads
`repeated_top_features__top30.csv`, selects the first 30 representative
features, applies median imputation, and then runs the FAMD -> KNN -> Louvain
-> UMAP workflow.

Large CV detail tables, raw feature libraries, cell-level outputs, and
patient-level feature matrices are excluded.
