# Limitations

This repository should be interpreted as a workflow and validation artifact, not as evidence of a clinically deployable predictor.

## Cohort Size

The analysis is based on 73 NSCLC patients treated with immune checkpoint inhibitors. This is a small sample size for a high-dimensional feature search over clinical and PBMC scRNA-seq-derived immune features.

## No External Validation Cohort

The retained results include internal 3-fold outer validation, but no independent external cohort. Without external validation, the findings should be considered exploratory.

## Protected Data Are Not Public

Raw clinical metadata, cell-level scRNA-seq outputs, and full patient-level feature matrices are excluded because of patient privacy and data-use restrictions. This limits full public reproducibility of the original analysis.

## Feature Search Multiplicity

The workflow scans many candidate immune features across multiple feature families (949 features). Even with nested cross-validation, correlation pruning, VIF pruning, and outer validation, the search space creates a substantial risk of false-positive or unstable feature selection.

## Fold Instability

Top-ranked features vary across outer folds. This instability is expected in a small cohort and should be interpreted as a warning against overclaiming any single biomarker or immune signature (see workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/outer_selected_candidates.csv, workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top30__details.csv).

## Negative Outer Validation Result

In the retained 3-fold outer validation, the single selected immune-feature model from each fold underperformed the clinical baseline model on held-out ROC-AUC and AUPRC. This does not mean every candidate feature was uninformative. Rather, it means the supervised feature-selection procedure did not produce a robust held-out predictive model in this small cohort (see workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/outer_fold_metrics.csv and workflows/results/feature_search_base_v2_single_feature_outer/outer_search_validation/outer_metrics_summary.csv).

## Exploratory Clustering Is Not Confirmatory

The FAMD -> KNN -> Louvain -> UMAP analysis using top-30 feature groups is exploratory. It is useful for visualizing possible patient subgroups and generating biological hypotheses, but it does not establish a validated patient stratification model.

## Intended Use

The repository is intended to demonstrate the design and implementation of a leakage-aware small-cohort translational ML workflow. It should not be used for clinical decision-making or treated as a validated predictive model.
