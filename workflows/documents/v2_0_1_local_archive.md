# Version 2.0.1 Local Archive Notes

Source folder:

```text
analysis/n73_manual_trial_20260305/version2/v.2.0.1
```

This local archive contains the v2.0.1 feature extraction scripts, feature
search scripts/config snapshots, and exploratory figures used after the
outer-validation feature search. The Git workflow keeps only reviewable scripts,
summary tables, and representative figures. Raw patient-level feature matrices
and cell-level objects are excluded by the repository data policy.

## Retained Summary Tables

The 3-fold outer-validation top-feature summaries are retained under:

```text
workflows/version2/results/feature_search/base_v2_single_feature_outer/outer_search_validation/
  repeated_top_features__top30.csv
  repeated_top_features__top30_details.csv
```

`repeated_top_features__top30.csv` summarizes the union of features appearing in
the top 30 within at least one of the 3 outer folds. The clustering branch uses
the first 30 representative features after sorting by mean group rank.

## Top-30 FAMD/KNN/Louvain/UMAP Branch

The exploratory clustering branch is implemented in:

```text
scripts/r/20260408_feature_based_clustering.R
```

Workflow:

1. Read `repeated_top_features__top30.csv`.
2. Select the first 30 representative features.
3. Merge those immune features with clinical covariates.
4. Median-impute missing numeric feature values.
5. Run FAMD on the mixed clinical + immune table.
6. Use FAMD dimensions 2-5 as the exploratory embedding.
7. Build a k-nearest-neighbor graph with `k = 5`.
8. Run Louvain clustering with `resolution = 1`.
9. Run UMAP with `n_neighbors = min(10, n - 1)`, `min_dist = 0.3`, and
   Euclidean distance.
10. Overlay cluster, binarized response, RECIST response, 6-month PFS, and
    PD-L1 TPS.

## Retained Figures

```text
workflows/version2/figures/exploratory_rplots/
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
  20260429_version2_trial1_FAMD_cluster_DEF.pdf
  20260429_version2_trial1_FAMD_cluster_DEF_categorical.pdf
  20260429_version2_trial1_FAMD_cluster_DEF_v2.pdf
```

## Excluded From Git

The following local archive content is intentionally not copied into Git:

- patient-level feature matrices such as `final_feature_df_with_cellrank_b3_filt.csv`
- raw or intermediate single-cell objects
- `.DS_Store`, `.af`, and other editor/native figure-editing sidecars
- large per-cell or per-patient intermediate tables
