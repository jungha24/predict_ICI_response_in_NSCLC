
# protocol note

1. Stage 1 feature library
    - Build the patient-level feature library first, then pass it into the Python search pipeline.
    - feature type:
        - cell-level composition
        - potency/dynamics or CellRank-like transition/priming surrogates
        - curated gene/pathway scores
        - pseudobulk de novo program/module
        - selected latent axes, e.g. PC summaries
    - major cell-type blocks:
        - B lineage
        - monocyte
        - NK
        - CD4 T
        - CD8 T
        - nonconventional T
    - QC/filtering:
        - retain only selected patients
        - apply minimum-patient, non-zero/detection, and unstable-column filters


    <details> <summary><b>Patient-level immune feature library</b></summary>
    
    | Level         | Feature category           | B cells                                                                                                             | Monocytes                                                                                           | NK cells                                                                                          | CD4 T cells                                                                                                        | CD8 T cells                                                                                                        | Non-conventional T cells                                                                                        |
    | ------------- | -------------------------- | ------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------- |
    | Cell level    | Composition                | Naive / memory                                                                                                      | -                                                                                                   | -                                                                                                 | resting naive / PDCD1hi activated / Treg                                                                           | Naive / cytotoxic                                                                                                  | -                                                                                                               |
    | Cell level    | PCA centroid-based balance | naive–memory                                                                                                        | active–non-classic                                                                                  | early response–resting                                                                            | resting naive–immediate early<br>resting naive–PDCD1hi activated<br>immediate early–PDCD1hi activated              | naive–cytotoxic                                                                                                    | TRDV2 TRGV9 gdT–TRDV1-like gdT                                                                                  |
    | Cell level    | Dynamics of cell status    | naive–memory–plasma                                                                                                 | active / non-classic                                                                                | immediate early / resting / CD56bright                                                            | resting naive / PDCD1hi activated / Treg                                                                           | naive / cytotoxic                                                                                                  | TRDV2 TRGV9 gdT / TRDV1-like gdT                                                                                |
    | Patient level | Curated gene sets          | 5 inflammation-related Hallmarks<br>3 antigen presentation-related gene sets<br>11 B cell-related functional panels | 5 inflammation-related Hallmarks<br>12 monocyte-related functional panels<br>4 additional gene sets | 5 inflammation-related Hallmarks<br>6 NK cell-related functional panels<br>4 additional gene sets | 10 inflammation- and status-related Hallmarks<br>14 CD4 T cell-related functional panels<br>9 additional gene sets | 10 inflammation- and status-related Hallmarks<br>14 CD8 T cell-related functional panels<br>9 additional gene sets | 10 inflammation- and status-related Hallmarks<br>14 T cell-related functional panels<br>10 additional gene sets |
    | Patient level | De novo gene sets          | 4 metaprograms                                                                                                      | 6 metaprograms                                                                                      | 3 metaprograms                                                                                    | 5 metaprograms                                                                                                     | 4 metaprograms                                                                                                     | 1 metaprogram                                                                                                   |
    
    </details>



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
src/feature_search_base_v2/io_utils.py # input output utils
src/feature_search_base_v2/search.py # run workflow
src/feature_search_base_v2/data.py # make clinical + immune feature matrix
src/feature_search_base_v2/design.py # define endpoint/baseline feature set
src/feature_search_base_v2/models.py # model engine
src/report_repeated_outer_top_features.py # extract top features per fold
```

search.py
  - Stage 2: run the baseline clinical model for each endpoint, add one candidate immune feature at a time, and save the delta metric versus baseline.
  - Stage 3: start from the top single-feature rankings, apply family caps, reduce redundant features with pairwise correlation and VIF pruning, and retain final candidates.

models.py
  - training-fold-only imputation, standardization, and one-hot encoding
  - zero-variance, optional corr/VIF pruning << pruning
  - nested CV elastic-net logistic/Cox
  - bootstrap stability


When outer validation is disabled, Stage 2/3 are run once on the full analysis cohort and produce CV-based feature rankings and pruned candidate lists.

When outer validation is enabled, the same Stage 2/3 search procedure is rerun independently inside each outer-training fold. A candidate feature is selected from the outer-training search only, then both the clinical baseline and clinical-plus-selected-feature models are fit on the outer-training patients and evaluated on the held-out outer-test patients. These held-out evaluations are written to outer_search_validation/outer_fold_metrics.csv. The repeated top-feature summaries aggregate feature-level ranking signals across the outer-training searches and are used for exploratory FAMD/KNN/Louvain/UMAP analysis, not as direct held-out performance metrics.

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
