# Reviewer Guide

This repository documents a small-cohort translational ML workflow for predicting immune checkpoint inhibitor response in NSCLC using clinical metadata and PBMC scRNA-seq-derived immune features.

## What to Review First

1. Pipeline code
   - `src/run_feature_search_v2.py`
   - `src/feature_search_base_v2/search.py`
   - `src/feature_search_base_v2/data.py`
   - `src/feature_search_base_v2/design.py`
   - `src/feature_search_base_v2/models.py`
   
2. Main configuration
   - `configs/version2_modeling_base_v2_single_feature_outer.yaml`
   - `workflows/configs/version2_modeling_base_v2_single_feature_outer.yaml`

3. Outer validation summary
   - `workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/outer_metrics_summary.csv`
   - `workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/outer_fold_metrics.csv`
   - `workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/outer_selected_candidates.csv`

4. Top feature stability across outer folds
   - `workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top30.csv`
   - `workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top30__details.csv`

5. Exploratory subgroup analysis
   - `scripts/r/20260408_feature_based_clustering.R`
   - `workflows/figures/exploratory_rplots/`

## Workflow Summary

The workflow starts from a patient-level immune feature library derived from PBMC scRNA-seq. Candidate immune features are tested as add-ons to a clinical baseline model. Features are ranked by nested cross-validation performance, then filtered using biological family caps, pairwise correlation pruning, and VIF-based redundancy pruning.

The retained version includes a 3-fold outer validation. In each outer fold, feature search is rerun only on the outer-training patients, and the selected candidate is evaluated on held-out outer-test patients.

The repository also includes exploratory FAMD -> KNN -> Louvain -> UMAP analysis using the top-30 feature groups summarized across the outer-training searches. This clustering analysis is intended for hypothesis generation and visualization, not confirmatory model validation.

## How to Interpret the Results

The feature search and the final selected-model validation should be read separately. 
- The top-30 feature table captures feature-level signal and stability across outer-training searches, and it is used for exploratory subgroup analysis.
- The final supervised outer-validation result evaluates the single selected immune-feature model from each fold on held-out patients; in this retained run, those selected models underperformed the clinical baseline on held-out ROC-AUC and AUPRC. This supports the conclusion that the cohort is too small and/or the feature space too large for a stable predictive claim without additional validation.

The main value of this repository is the workflow structure: protected clinical data handling, staged feature search, leakage-aware validation, redundancy pruning, retained run artifacts, and explicit documentation of negative validation evidence.
