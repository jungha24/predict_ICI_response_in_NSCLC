
# protocol note

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
   - Output: nested CV delta metric per endpoint
        - inner CV: modeling.resampling.n_inner_splits: 4
        - random seed: modeling.random_state: 42
        - logistic elastic-net grid: modeling.logistic.alpha_grid, l1_ratio_grid
        - ranking metric: search.stage2.ranking_metric_binary: delta_roc_auc_mean
        - baseline: search.baseline.mode: clinical_with_pd_l1
3. Stage 3: redundancy pruning
   - Biological family caps
   - Pairwise correlation pruning
   - VIF pruning

4. Outer validation
   - The run attempted one endpoint, `Binarized_response`, and generated 3 outer folds in total. In each outer fold, the feature search was rerun on the outer-training patients, then the selected candidate was evaluated on the held-out outer-test patients.
        - `validation.outer_search_cv.enabled: true`
        - `n_splits: 3`
        - `n_repeats: 1`

    - Outer validation outputs:
        - `outer_search_validation/outer_validation_manifest.json`: executed settings and total fold count
        - `outer_search_validation/outer_selected_candidates.csv`: selected candidate feature per outer fold
        - `outer_search_validation/outer_fold_metrics.csv`: baseline and selected-model performance on each outer test fold
        - `outer_search_validation/binarized_response/fold_*/stage2_single_feature_scan.csv`: fold-specific stage2 feature ranking

5. Exploratory patient clustering from outer-validation features
   - Summarize the top 30 feature groups from each of the 3 outer-validation
     folds.
   - Use those features to build a patient-level exploratory matrix (Median-impute missing numeric feature values).
   - Run FAMD on the mixed clinical + immune feature matrix.
   - Build a k-nearest-neighbor graph on selected FAMD dimensions.
   - Run Louvain clustering and visualize the result with UMAP.
   - Overlay cluster, binarized response, RECIST response, 6-month PFS, and
     PD-L1 TPS on the UMAP.

## Main Code

```text
src/run_feature_search_v2.py #entrypoint
src/feature_search_base_v2/io_utils.py #input output utils
src/feature_search_base_v2/data.py #integrate clinical/feature matrix
src/feature_search_base_v2/design.py #define design
src/feature_search_base_v2/models.py #model engine
src/feature_search_base_v2/search.py #실제 탐색 흐름
src/report_repeated_outer_top_features.py # extract top features per fold
```

models.py
  - training-fold 기준 imputation/standardization/one-hot
  - zero-variance, optional corr/VIF pruning << pruning
  - nested CV elastic-net logistic/Cox
  - bootstrap stability

search.py
  - stage 2: baseline clinical model을 endpoint별로 돌리고, candidate immune feature을 하나씩 더한 걸 전부 평가해 baseline 대비 delta metric저장
  - stage 3: single-feature ranking 상위권에서 family cap을 걸고, pairwise corr와 VIF로 중복 feature 줄여 최종 후보 남기기 << pruning
  
## Key Configs

```text
configs/version2_modeling_smoketest.yaml
configs/version2_modeling_base.yaml
configs/version2_modeling_base_v2.yaml
configs/version2_modeling_base_v2_single_feature_outer.yaml
```

`version2_modeling_base_v2_single_feature_outer.yaml` is the main retained

## Representative Command

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src \
conda run -n nsclc-subtype python src/version2/run_feature_search_v2.py \
  --config configs/version2_modeling_base_v2_single_feature_outer.yaml

conda run -n nsclc-subtype python src/version2/report_repeated_outer_top_features.py \
  --run-dir data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer \
  --top-n 30
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
workflows/results/
  stage2_baseline_metrics.csv
  stage2_single_feature_scan.csv
  stage3_pruned_candidates.csv
  stage3_redundancy_pairs.csv
  stage3_vif_trace.csv
  outer_search_validation/
    repeated_top_features__top30.csv
    repeated_top_features__top30_details.csv
```

The clustering script is retained at
`scripts/r/20260408_feature_based_clustering.R`. Its top-30 branch reads
`repeated_top_features__top30.csv`, selects the first 30 representative
features, applies median imputation, and then runs the FAMD -> KNN -> Louvain
-> UMAP workflow.

Cell-level outputs, and patient-level feature matrices are excluded.
