#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Allow direct execution via `python src/version1/run_subtype_pipeline.py`.
SRC_ROOT = Path(__file__).resolve().parents[1]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from version1.subtype.cohort import build_discovery_cohort, summarize_discovery_cohort
from version1.subtype.features import (
    build_patient_metadata,
    fill_missing_values,
    make_proportion_features,
)
from version1.subtype.pseudobulk import aggregate_pseudobulk, get_valid_pseudobulk_groups
from version1.subtype.feature_selection import (
    check_feature_residual_signal,
    check_pca_nmf_redundancy,
    check_prop_related_features,
    select_features_for_clustering,
)
from version1.subtype.external_features import load_external_patient_features
from version1.subtype.batch_adjust import resolve_covariates
from version1.subtype.clustering import run_patient_clustering
from version1.subtype.io_utils import ensure_outdir, load_inputs, read_yaml, write_json, write_yaml
from version1.subtype.plotting import save_basic_plots
from version1.subtype.scanpy_pseudobulk import run_combat_seq_nmf_workflow, run_scanpy_pca_workflow


def log_step(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def collapse_group_level_nmf_to_patient(nmf_module_scores_df: pd.DataFrame | None) -> pd.DataFrame | None:
    # Convert group-level NMF module tables back to one patient-level row for diagnostics.
    if nmf_module_scores_df is None or nmf_module_scores_df.empty:
        return None
    value_columns = [col for col in nmf_module_scores_df.columns if col not in {"analysis_id", "group_key"}]
    if not value_columns:
        return None
    return nmf_module_scores_df.groupby("analysis_id", as_index=False)[value_columns].first()


def main(config_path: str):
    log_step(f"Starting subtype pipeline with config: {config_path}")
    cfg = read_yaml(config_path)
    outdir = ensure_outdir(cfg["paths"]["outdir"])

    log_step(f"Preparing output directory: {outdir}")
    write_yaml(cfg, outdir / "config.yaml")

    log_step("Loading input matrices and metadata")
    counts, genes, barcodes, meta = load_inputs(cfg)

    log_step("Building discovery cohort")
    counts_use, meta_use, cohort_info = build_discovery_cohort(counts, meta, cfg)
    cohort_summary = summarize_discovery_cohort(meta_use)

    log_step("Building patient-level metadata plus configured immune feature inputs")
    patient_meta = build_patient_metadata(meta_use)
    prop_wide, prop_long, prop_diag_wide = make_proportion_features(meta_use, cfg)
    external_feature_df, external_feature_manifest = load_external_patient_features(patient_meta, cfg)

    log_step("Selecting valid pseudobulk groups and aggregating counts")
    valid_groups = get_valid_pseudobulk_groups(meta_use, cfg)
    pb_counts, pb_group_names, pb_meta = aggregate_pseudobulk(counts_use, meta_use, valid_groups)
    if valid_groups.empty:
        log_step("No valid pseudobulk groups found; pseudobulk-based PCA/NMF features will be skipped")

    batch_cfg = cfg.get("batch_adjustment", {})
    pca_cfg = cfg.get("pca_features", {})
    covariates = resolve_covariates(cfg) if batch_cfg.get("enabled", True) else []
    pca_gene_set = "Scanpy HVGs" if pca_cfg.get("use_hvg", True) else "all nonzero-variance genes"
    if covariates:
        log_step(f"Running Scanpy pseudobulk PCA on {pca_gene_set} with gene-level corrected logCPM")
    else:
        log_step(f"Running Scanpy pseudobulk PCA on {pca_gene_set} without batch correction")
    pca_feature_df, pca_variance_df = run_scanpy_pca_workflow(
        pb_counts,
        pb_group_names,
        pb_meta,
        genes,
        cfg,
        outdir=outdir,
        covariates=covariates,
    )

    log_step("Running ComBat-seq NMF review workflow")
    nmf_feature_df, nmf_module_scores_df, nmf_meta, nmf_gene_list_df = run_combat_seq_nmf_workflow(
        pb_counts,
        pb_group_names,
        pb_meta,
        genes,
        cfg,
        outdir=outdir,
        covariates=covariates,
    )

    log_step("Merging patient feature table and filling missing values")
    nmf_diagnostic_feature_df = collapse_group_level_nmf_to_patient(nmf_module_scores_df)
    patient_features = patient_meta.merge(prop_wide, on="analysis_id", how="left")
    diagnostic_features = patient_meta.merge(prop_diag_wide, on="analysis_id", how="left")
    if external_feature_df.shape[1] > 1:
        patient_features = patient_features.merge(external_feature_df, on="analysis_id", how="left")
        diagnostic_features = diagnostic_features.merge(external_feature_df, on="analysis_id", how="left")
    if pca_feature_df is not None:
        patient_features = patient_features.merge(pca_feature_df, on="analysis_id", how="left")
        diagnostic_features = diagnostic_features.merge(pca_feature_df, on="analysis_id", how="left")
    if nmf_feature_df is not None:
        patient_features = patient_features.merge(nmf_feature_df, on="analysis_id", how="left")
        diagnostic_features = diagnostic_features.merge(nmf_feature_df, on="analysis_id", how="left")
    if nmf_diagnostic_feature_df is not None:
        diagnostic_features = diagnostic_features.merge(nmf_diagnostic_feature_df, on="analysis_id", how="left")

    patient_features = fill_missing_values(patient_features)
    diagnostic_features = fill_missing_values(diagnostic_features)

    log_step("Running feature diagnostics and selecting clustering features")
    prop_related_df = check_prop_related_features(diagnostic_features, cfg)
    resid_signal_df = check_feature_residual_signal(diagnostic_features, cfg)
    pca_nmf_red_df = check_pca_nmf_redundancy(diagnostic_features)

    selected_features, dropped_features = select_features_for_clustering(
        patient_features,
        resid_signal_df,
        pca_nmf_red_df,
        cfg,
    )

    patient_features_adj = patient_features.copy()

    log_step("Running patient clustering and dimensionality reduction")
    pca_df, umap_df, cluster_df, sil_df, metrics, patient_pca_variance_df = run_patient_clustering(
        patient_features_adj,
        selected_features,
        cfg,
    )

    patient_final = patient_features_adj.merge(cluster_df, on="analysis_id", how="left")

    log_step("Saving plots")
    save_basic_plots(outdir, pca_df, umap_df, patient_final, patient_pca_variance_df, cfg)

    log_step("Writing output tables and metadata files")
    patient_features.to_csv(outdir / "patient_features.csv", index=False)
    external_patient_features_out = patient_meta[["analysis_id", "patient_id"]].merge(
        external_feature_df, on="analysis_id", how="left"
    )
    external_patient_features_out.to_csv(outdir / "external_patient_features.csv", index=False)
    prop_long.to_csv(outdir / "celltype_proportions_long.csv", index=False)
    valid_groups.to_csv(outdir / "valid_pseudobulk_groups.csv", index=False)
    pb_meta.to_csv(outdir / "pseudobulk_group_metadata.csv", index=False)
    pca_variance_df.to_csv(outdir / "pca_explained_variance.csv", index=False)
    if nmf_module_scores_df is not None:
        nmf_module_scores_df.to_csv(outdir / "nmf_module_scores.csv", index=False)
    else:
        pd.DataFrame(columns=["analysis_id", "group_key"]).to_csv(outdir / "nmf_module_scores.csv", index=False)
    nmf_gene_list_df.to_csv(outdir / "nmf_gene_lists.csv", index=False)
    prop_related_df.to_csv(outdir / "prop_related_df.csv", index=False)
    resid_signal_df.to_csv(outdir / "resid_signal_df.csv", index=False)
    pca_nmf_red_df.to_csv(outdir / "pca_nmf_redundancy.csv", index=False)
    pca_df.to_csv(outdir / "patient_pca.csv", index=False)
    patient_pca_variance_df.to_csv(outdir / "patient_pca_explained_variance.csv", index=False)
    umap_df.to_csv(outdir / "patient_umap.csv", index=False)
    cluster_df.to_csv(outdir / "patient_clusters.csv", index=False)
    sil_df.to_csv(outdir / "silhouette_by_k.csv", index=False)

    write_json(
        {
            "selected_features": selected_features,
            "dropped_features": dropped_features,
            "covariates_used_for_adjustment": covariates,
            "batch_adjustment_stage": "pseudobulk_gene_level" if covariates else "disabled",
            "proportion_transform": (
                cfg.get("proportion_features", {}).get("transform", "clr")
                if cfg.get("proportion_features", {}).get("enabled", True)
                else "disabled"
            ),
            "external_patient_features_enabled": bool(external_feature_manifest.get("enabled", False)),
            "external_patient_features_path": external_feature_manifest.get("path"),
            "external_patient_feature_columns": external_feature_manifest.get("feature_columns", []),
            "pca_gene_set": "scanpy_hvg" if pca_cfg.get("use_hvg", True) else "all_nonzero_variance_genes",
            "pca_hvg_n_top_genes": pca_cfg.get("hvg_n_top_genes", 2000),
            "pca_hvg_flavor": pca_cfg.get("hvg_flavor", "seurat"),
            "patient_n_pcs_requested": cfg.get("patient_clustering", {}).get("n_pcs_patient"),
        },
        outdir / "selected_features.json",
    )
    write_json(nmf_meta, outdir / "nmf_top_genes.json")
    write_json(cohort_info, outdir / "cohort_info.json")
    write_json(cohort_summary, outdir / "cohort_summary.json")
    write_json(metrics, outdir / "metrics.json")

    log_step("Subtype pipeline completed successfully")
    print("Done.")
    print(json.dumps(metrics, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to subtype YAML config")
    args = parser.parse_args()
    main(args.config)
