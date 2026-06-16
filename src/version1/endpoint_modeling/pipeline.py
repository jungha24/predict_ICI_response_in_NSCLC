"""High-level orchestration for endpoint modeling analyses and output writing."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import helmert

from .data import build_prepared_dataset
from .design import AnalysisSpec, build_analysis_specs, resolve_endpoint_specs
from .io_utils import ensure_outdir, log_step, read_yaml, write_json, write_yaml


warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


def resolve_analysis_columns(
    spec: AnalysisSpec,
    clinical_feature_sets: dict[str, list[str]],
    immune_blocks: dict[str, object],
) -> tuple[list[str], list[str]]:
    # Each analysis mixes one named clinical set with zero or more immune blocks.
    clinical_cols = clinical_feature_sets.get(spec.clinical_feature_set, []) if spec.clinical_feature_set else []
    immune_cols: list[str] = []
    for block_name in spec.immune_blocks:
        if block_name not in immune_blocks:
            raise ValueError(f"Analysis '{spec.name}' references unknown immune block '{block_name}'.")
        immune_cols.extend(immune_blocks[block_name].columns)
    all_columns = list(dict.fromkeys(clinical_cols + immune_cols))
    return all_columns, clinical_cols


def build_execution_plan(
    analysis_specs_by_block: dict[str, list[AnalysisSpec]],
    clinical_feature_sets: dict[str, list[str]],
    immune_blocks: dict[str, object],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    # Identical analysis designs across primary/secondary blocks should be fit once and reused.
    execution_by_signature: dict[tuple[tuple[str, ...], tuple[str, ...]], dict[str, object]] = {}
    analysis_registry_rows: list[dict[str, object]] = []

    for block_name, specs in analysis_specs_by_block.items():
        for spec in specs:
            feature_cols, clinical_cols = resolve_analysis_columns(
                spec=spec,
                clinical_feature_sets=clinical_feature_sets,
                immune_blocks=immune_blocks,
            )
            if not feature_cols:
                analysis_registry_rows.append(
                    {
                        "analysis_block": block_name,
                        "feature_set": spec.name,
                        "clinical_feature_set": spec.clinical_feature_set,
                        "immune_blocks": ",".join(spec.immune_blocks),
                        "n_clinical_features": len(clinical_cols),
                        "n_immune_features": 0,
                        "description": spec.description,
                        "execution_block": pd.NA,
                        "execution_feature_set": pd.NA,
                        "is_execution_alias": pd.NA,
                    }
                )
                continue

            signature = (tuple(clinical_cols), tuple(feature_cols))
            if signature not in execution_by_signature:
                execution_by_signature[signature] = {
                    "canonical_block": block_name,
                    "canonical_spec": spec,
                    "feature_cols": feature_cols,
                    "clinical_cols": clinical_cols,
                    "aliases": [],
                }
            execution_by_signature[signature]["aliases"].append((block_name, spec))
            canonical_block = str(execution_by_signature[signature]["canonical_block"])
            canonical_spec = execution_by_signature[signature]["canonical_spec"]

            analysis_registry_rows.append(
                {
                    "analysis_block": block_name,
                    "feature_set": spec.name,
                    "clinical_feature_set": spec.clinical_feature_set,
                    "immune_blocks": ",".join(spec.immune_blocks),
                    "n_clinical_features": len(clinical_cols),
                    "n_immune_features": len([col for col in feature_cols if col not in clinical_cols]),
                    "description": spec.description,
                    "execution_block": canonical_block,
                    "execution_feature_set": canonical_spec.name,
                    "is_execution_alias": not (canonical_block == block_name and canonical_spec.name == spec.name),
                }
            )

    return list(execution_by_signature.values()), analysis_registry_rows


def build_comparison_summary(metrics_df: pd.DataFrame, comparisons_cfg: list[dict]) -> pd.DataFrame:
    # Post-hoc metric deltas summarize incremental value comparisons configured in YAML.
    if metrics_df.empty or not comparisons_cfg:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    for comparison in comparisons_cfg:
        block_name = comparison.get("analysis_block")
        reference = comparison["reference"]
        candidate = comparison["candidate"]

        ref_df = metrics_df[(metrics_df["analysis_block"] == block_name) & (metrics_df["feature_set"] == reference)]
        cand_df = metrics_df[(metrics_df["analysis_block"] == block_name) & (metrics_df["feature_set"] == candidate)]
        if ref_df.empty or cand_df.empty:
            continue

        merged = cand_df.merge(
            ref_df,
            on=["analysis_block", "endpoint"],
            suffixes=("_candidate", "_reference"),
        )
        numeric_columns = [
            column.replace("_candidate", "")
            for column in merged.columns
            if column.endswith("_candidate") and pd.api.types.is_numeric_dtype(merged[column])
        ]
        for row in merged.itertuples(index=False):
            record = {
                "comparison_name": comparison.get("name", f"{candidate}_vs_{reference}"),
                "analysis_block": getattr(row, "analysis_block"),
                "endpoint": getattr(row, "endpoint"),
                "candidate": candidate,
                "reference": reference,
            }
            for base_name in numeric_columns:
                candidate_value = getattr(row, f"{base_name}_candidate")
                reference_value = getattr(row, f"{base_name}_reference")
                record[f"{base_name}_candidate"] = candidate_value
                record[f"{base_name}_reference"] = reference_value
                if pd.notna(candidate_value) and pd.notna(reference_value):
                    record[f"delta_{base_name}"] = candidate_value - reference_value
                else:
                    record[f"delta_{base_name}"] = pd.NA
            rows.append(record)

    return pd.DataFrame(rows)


def build_ilr_block_metadata(prepared) -> dict[str, dict[str, object]]:
    metadata: dict[str, dict[str, object]] = {}
    feature_dictionary = prepared.feature_dictionary.copy()
    for block_name, block in prepared.immune_blocks.items():
        if str(block.transform).lower() != "ilr":
            continue
        block_rows = feature_dictionary[
            (feature_dictionary["block"] == block_name) & (feature_dictionary["feature_column"].isin(block.columns))
        ].copy()
        if block_rows.empty:
            continue
        source_columns_value = (
            block_rows["source_columns"].dropna().astype(str).iloc[0]
            if block_rows["source_columns"].dropna().shape[0] > 0
            else ""
        )
        part_order = [value.strip() for value in source_columns_value.split("|") if value.strip()]
        if len(part_order) < 2:
            continue
        metadata[block_name] = {
            "raw_features": list(block.columns),
            "part_order": part_order,
            "prefix": block.prefix,
        }
    return metadata


def build_ilr_celltype_weight_summary(
    coefficients_df: pd.DataFrame,
    ilr_block_metadata: dict[str, dict[str, object]],
) -> pd.DataFrame:
    if coefficients_df.empty or not ilr_block_metadata:
        return pd.DataFrame(
            columns=[
                "analysis_block",
                "endpoint",
                "feature_set",
                "ilr_block",
                "celltype",
                "clr_weight",
                "abs_clr_weight",
                "n_parts",
                "n_ilr_axes",
                "n_present_ilr_axes",
            ]
        )

    rows: list[dict[str, object]] = []
    group_cols = ["analysis_block", "endpoint", "feature_set"]
    for keys, group_df in coefficients_df.groupby(group_cols, dropna=False):
        analysis_block, endpoint, feature_set = keys
        raw_coef_map = group_df.groupby("raw_feature", dropna=False)["coef"].first().to_dict()
        for block_name, metadata in ilr_block_metadata.items():
            raw_features = list(metadata["raw_features"])
            part_order = list(metadata["part_order"])
            basis = helmert(len(part_order), full=False).T
            beta_ilr = np.zeros(len(raw_features), dtype=float)
            n_present = 0
            for idx, raw_feature in enumerate(raw_features):
                if raw_feature in raw_coef_map and pd.notna(raw_coef_map[raw_feature]):
                    beta_ilr[idx] = float(raw_coef_map[raw_feature])
                    n_present += 1
            if n_present == 0:
                continue

            clr_weights = basis @ beta_ilr
            for celltype, clr_weight in zip(part_order, clr_weights.tolist()):
                rows.append(
                    {
                        "analysis_block": analysis_block,
                        "endpoint": endpoint,
                        "feature_set": feature_set,
                        "ilr_block": block_name,
                        "celltype": celltype,
                        "clr_weight": float(clr_weight),
                        "abs_clr_weight": float(abs(clr_weight)),
                        "n_parts": len(part_order),
                        "n_ilr_axes": len(raw_features),
                        "n_present_ilr_axes": n_present,
                    }
                )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(
        ["analysis_block", "endpoint", "feature_set", "ilr_block", "abs_clr_weight"],
        ascending=[True, True, True, True, False],
    ).reset_index(drop=True)


def run_endpoint_modeling(config_path: str) -> Path:
    try:
        # Delay heavy imports so config/data utilities still work even if modeling deps are missing.
        from .models import (
            bootstrap_stability_cox,
            bootstrap_stability_logistic,
            evaluate_cox_nested_cv,
            evaluate_logistic_nested_cv,
            tune_on_full_data_cox,
            tune_on_full_data_logistic,
        )
    except ModuleNotFoundError as exc:
        missing_name = getattr(exc, "name", None) or str(exc)
        raise ModuleNotFoundError(
            "Endpoint modeling requires optional packages including `scikit-learn` and `lifelines`. "
            f"Missing dependency: {missing_name}. Install or activate the intended modeling environment first."
        ) from exc

    cfg = read_yaml(config_path)
    outdir = ensure_outdir(cfg["paths"]["output_dir"])

    log_step(f"Starting endpoint modeling with config: {config_path}")
    write_yaml(cfg, outdir / "config.yaml")

    # Step 1: build one matched analysis table and record which source features survived QC.
    log_step(
        "Loading clinical metadata and patient-level immune features "
        f"(clinical={cfg['clinical']['path']}, "
        f"immune_feature_table={cfg['immune_features'].get('feature_table_path', cfg['immune_features'].get('path'))})"
    )
    prepared = build_prepared_dataset(cfg)
    prepared.analysis_df.to_csv(outdir / "merged_analysis_dataset.tsv", sep="\t", index=False)
    prepared.feature_dictionary.to_csv(outdir / "feature_dictionary.csv", index=False)
    for name, table in prepared.diagnostic_tables.items():
        table.to_csv(outdir / f"{name}.csv", index=False)
    write_json(prepared.manifest, outdir / "input_manifest.json")

    # Step 2: resolve the primary/secondary model plan from config.
    analysis_specs_by_block = build_analysis_specs(cfg)
    binary_endpoints, survival_endpoints = resolve_endpoint_specs(cfg, prepared.analysis_df.columns.tolist())

    endpoint_registry_rows: list[dict[str, object]] = []
    for endpoint in binary_endpoints:
        endpoint_df = prepared.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
        endpoint_registry_rows.append(
            {
                "endpoint_name": endpoint.name,
                "endpoint_type": "binary",
                "outcome_col": endpoint.outcome_col,
                "n_patients": int(endpoint_df.shape[0]),
                "n_events": int(endpoint_df[endpoint.outcome_col].sum()) if not endpoint_df.empty else 0,
                "n_nonevents": int((1 - endpoint_df[endpoint.outcome_col]).sum()) if not endpoint_df.empty else 0,
                "description": endpoint.description,
            }
        )
    for endpoint in survival_endpoints:
        endpoint_df = prepared.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
        endpoint_registry_rows.append(
            {
                "endpoint_name": endpoint.name,
                "endpoint_type": "survival",
                "time_col": endpoint.time_col,
                "event_col": endpoint.event_col,
                "n_patients": int(endpoint_df.shape[0]),
                "n_events": int(endpoint_df[endpoint.event_col].sum()) if not endpoint_df.empty else 0,
                "description": endpoint.description,
            }
        )

    resampling_cfg = cfg.get("modeling", {}).get("resampling", {})
    stability_cfg = cfg.get("modeling", {}).get("stability", {})
    logistic_cfg = cfg.get("modeling", {}).get("logistic", {})
    cox_cfg = cfg.get("modeling", {}).get("cox", {})
    collinearity_cfg = cfg.get("modeling", {}).get("multicollinearity", {})

    random_state = int(cfg.get("modeling", {}).get("random_state", 42))
    n_outer_splits = int(resampling_cfg.get("n_outer_splits", 5))
    n_outer_repeats = int(resampling_cfg.get("n_outer_repeats", 10))
    n_inner_splits = int(resampling_cfg.get("n_inner_splits", 4))
    n_bootstrap = int(stability_cfg.get("n_bootstrap", 200))
    selection_threshold = float(stability_cfg.get("selection_threshold", 1e-8))
    selection_scope = str(stability_cfg.get("selection_scope", "immune_only"))

    logistic_alpha_grid = list(logistic_cfg.get("alpha_grid", [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]))
    logistic_l1_grid = list(logistic_cfg.get("l1_ratio_grid", [0.0, 0.25, 0.5, 0.75, 1.0]))
    logistic_max_iter = int(logistic_cfg.get("max_iter", 5000))
    logistic_selection_rule = str(logistic_cfg.get("selection_rule", "best"))

    cox_penalizer_grid = list(cox_cfg.get("penalizer_grid", [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]))
    cox_l1_grid = list(cox_cfg.get("l1_ratio_grid", [0.0, 0.25, 0.5, 0.75, 1.0]))
    cox_clinical_penalty_factor = float(cox_cfg.get("clinical_penalty_factor", 0.01))
    cox_selection_rule = str(cox_cfg.get("selection_rule", "best"))

    execution_plan, analysis_registry_rows = build_execution_plan(
        analysis_specs_by_block=analysis_specs_by_block,
        clinical_feature_sets=prepared.clinical_feature_sets,
        immune_blocks=prepared.immune_blocks,
    )

    metrics_rows: list[dict] = []
    tuning_rows: list[dict] = []
    coefficient_frames: list[pd.DataFrame] = []

    # Step 3: run every planned analysis for all configured endpoints.
    for execution in execution_plan:
        block_name = str(execution["canonical_block"])
        spec = execution["canonical_spec"]
        feature_cols = list(execution["feature_cols"])
        clinical_cols = list(execution["clinical_cols"])
        aliases = list(execution["aliases"])
        if not feature_cols:
            log_step(f"Skipping {block_name}/{spec.name} because no usable features remain after QC")
            continue

        if len(aliases) > 1:
            alias_labels = ", ".join(f"{alias_block}/{alias_spec.name}" for alias_block, alias_spec in aliases[1:])
            log_step(f"Reusing one execution for duplicate analysis designs: {block_name}/{spec.name} -> {alias_labels}")

        # Step 3a: binary endpoints -> repeated-CV elastic-net logistic + bootstrap stability.
        for endpoint in binary_endpoints:
            endpoint_df = prepared.analysis_df.dropna(subset=[endpoint.outcome_col]).copy()
            if endpoint_df.empty or endpoint_df[endpoint.outcome_col].nunique(dropna=True) < 2:
                log_step(
                    f"Skipping {block_name}/{spec.name} for {endpoint.name} because fewer than two outcome classes are available"
                )
                continue

            log_step(f"Running {block_name}/{spec.name} for {endpoint.name}")
            try:
                log_fold_df, log_summary = evaluate_logistic_nested_cv(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    outcome_col=endpoint.outcome_col,
                    alpha_grid=logistic_alpha_grid,
                    l1_ratio_grid=logistic_l1_grid,
                    n_outer_splits=n_outer_splits,
                    n_outer_repeats=n_outer_repeats,
                    n_inner_splits=n_inner_splits,
                    random_state=random_state,
                    max_iter=logistic_max_iter,
                    selection_rule=logistic_selection_rule,
                    collinearity_cfg=collinearity_cfg,
                )
                log_alpha, log_l1_ratio, log_inner_auc, log_coef_df = tune_on_full_data_logistic(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    outcome_col=endpoint.outcome_col,
                    alpha_grid=logistic_alpha_grid,
                    l1_ratio_grid=logistic_l1_grid,
                    n_inner_splits=n_inner_splits,
                    random_state=random_state,
                    max_iter=logistic_max_iter,
                    selection_rule=logistic_selection_rule,
                    collinearity_cfg=collinearity_cfg,
                )
                log_stability = bootstrap_stability_logistic(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    outcome_col=endpoint.outcome_col,
                    alpha=log_alpha,
                    l1_ratio=log_l1_ratio,
                    n_bootstrap=n_bootstrap,
                    selection_threshold=selection_threshold,
                    random_state=random_state,
                    max_iter=logistic_max_iter,
                    selection_scope=selection_scope,
                    collinearity_cfg=collinearity_cfg,
                )
            except ValueError as exc:
                log_step(f"Skipping {block_name}/{spec.name} for {endpoint.name}: {exc}")
                continue

            for alias_block, alias_spec in aliases:
                log_fold_df.to_csv(outdir / f"{alias_block}__{alias_spec.name}__{endpoint.slug}_cv_folds.csv", index=False)
                log_stability.to_csv(outdir / f"{alias_block}__{alias_spec.name}__{endpoint.slug}_stability.csv", index=False)
                metrics_rows.append(
                    {
                        "analysis_block": alias_block,
                        "endpoint": endpoint.name,
                        "feature_set": alias_spec.name,
                        "n_endpoint_patients": int(endpoint_df.shape[0]),
                        **log_summary,
                    }
                )
                tuning_rows.append(
                    {
                        "analysis_block": alias_block,
                        "endpoint": endpoint.name,
                        "feature_set": alias_spec.name,
                        "best_alpha": log_alpha,
                        "best_l1_ratio": log_l1_ratio,
                        "full_inner_auc": log_inner_auc,
                        "n_features_after_preprocessing": int(log_coef_df.shape[0]),
                        "selection_rule": logistic_selection_rule,
                        "logistic_force_in_note": str(logistic_cfg.get("force_in_note", "")),
                    }
                )
                alias_coef_df = log_coef_df.copy()
                alias_coef_df.insert(0, "feature_set", alias_spec.name)
                alias_coef_df.insert(0, "endpoint", endpoint.name)
                alias_coef_df.insert(0, "analysis_block", alias_block)
                coefficient_frames.append(alias_coef_df)

        # Step 3b: survival endpoints -> repeated-CV elastic-net Cox + bootstrap stability.
        for endpoint in survival_endpoints:
            endpoint_df = prepared.analysis_df.dropna(subset=[endpoint.time_col, endpoint.event_col]).copy()
            if endpoint_df.empty or endpoint_df[endpoint.event_col].nunique(dropna=True) < 2:
                log_step(
                    f"Skipping {block_name}/{spec.name} for {endpoint.name} because fewer than two event classes are available"
                )
                continue

            log_step(f"Running {block_name}/{spec.name} for {endpoint.name}")
            try:
                cox_fold_df, cox_summary = evaluate_cox_nested_cv(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    time_col=endpoint.time_col,
                    event_col=endpoint.event_col,
                    penalizer_grid=cox_penalizer_grid,
                    l1_ratio_grid=cox_l1_grid,
                    n_outer_splits=n_outer_splits,
                    n_outer_repeats=n_outer_repeats,
                    n_inner_splits=n_inner_splits,
                    random_state=random_state,
                    clinical_penalty_factor=cox_clinical_penalty_factor,
                    selection_rule=cox_selection_rule,
                    collinearity_cfg=collinearity_cfg,
                )
                cox_penalizer, cox_l1_ratio, cox_inner_cindex, cox_coef_df = tune_on_full_data_cox(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    time_col=endpoint.time_col,
                    event_col=endpoint.event_col,
                    penalizer_grid=cox_penalizer_grid,
                    l1_ratio_grid=cox_l1_grid,
                    n_inner_splits=n_inner_splits,
                    random_state=random_state,
                    clinical_penalty_factor=cox_clinical_penalty_factor,
                    selection_rule=cox_selection_rule,
                    collinearity_cfg=collinearity_cfg,
                )
                cox_stability = bootstrap_stability_cox(
                    df=endpoint_df,
                    feature_cols=feature_cols,
                    clinical_cols=clinical_cols,
                    time_col=endpoint.time_col,
                    event_col=endpoint.event_col,
                    penalizer=cox_penalizer,
                    l1_ratio=cox_l1_ratio,
                    n_bootstrap=n_bootstrap,
                    selection_threshold=selection_threshold,
                    random_state=random_state,
                    clinical_penalty_factor=cox_clinical_penalty_factor,
                    selection_scope=selection_scope,
                    collinearity_cfg=collinearity_cfg,
                )
            except ValueError as exc:
                log_step(f"Skipping {block_name}/{spec.name} for {endpoint.name}: {exc}")
                continue

            for alias_block, alias_spec in aliases:
                cox_fold_df.to_csv(outdir / f"{alias_block}__{alias_spec.name}__{endpoint.slug}_cv_folds.csv", index=False)
                cox_stability.to_csv(outdir / f"{alias_block}__{alias_spec.name}__{endpoint.slug}_stability.csv", index=False)
                metrics_rows.append(
                    {
                        "analysis_block": alias_block,
                        "endpoint": endpoint.name,
                        "feature_set": alias_spec.name,
                        "n_endpoint_patients": int(endpoint_df.shape[0]),
                        **cox_summary,
                    }
                )
                tuning_rows.append(
                    {
                        "analysis_block": alias_block,
                        "endpoint": endpoint.name,
                        "feature_set": alias_spec.name,
                        "best_penalizer": cox_penalizer,
                        "best_l1_ratio": cox_l1_ratio,
                        "full_inner_cindex": cox_inner_cindex,
                        "n_features_after_preprocessing": int(cox_coef_df.shape[0]),
                        "selection_rule": cox_selection_rule,
                        "cox_clinical_penalty_factor": cox_clinical_penalty_factor,
                    }
                )
                alias_coef_df = cox_coef_df.copy()
                alias_coef_df.insert(0, "feature_set", alias_spec.name)
                alias_coef_df.insert(0, "endpoint", endpoint.name)
                alias_coef_df.insert(0, "analysis_block", alias_block)
                coefficient_frames.append(alias_coef_df)

    # Step 4: consolidate summary tables across all analyses/endpoints.
    endpoint_registry_df = pd.DataFrame(endpoint_registry_rows)
    endpoint_registry_df.to_csv(outdir / "endpoint_registry.csv", index=False)

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(outdir / "model_metrics_summary.csv", index=False)

    tuning_df = pd.DataFrame(tuning_rows)
    tuning_df.to_csv(outdir / "full_data_tuning_summary.csv", index=False)

    analysis_registry_df = pd.DataFrame(analysis_registry_rows)
    analysis_registry_df.to_csv(outdir / "analysis_registry.csv", index=False)

    coefficients_df = pd.concat(coefficient_frames, ignore_index=True) if coefficient_frames else pd.DataFrame()
    if not coefficients_df.empty:
        coefficients_df = coefficients_df.sort_values(
            ["analysis_block", "endpoint", "feature_set", "feature_scope", "abs_coef"],
            ascending=[True, True, True, True, False],
        )
    coefficients_df.to_csv(outdir / "full_data_coefficients.csv", index=False)

    ilr_block_metadata = build_ilr_block_metadata(prepared)
    ilr_weight_df = build_ilr_celltype_weight_summary(coefficients_df, ilr_block_metadata)
    ilr_weight_df.to_csv(outdir / "full_data_ilr_celltype_weights.csv", index=False)

    comparison_df = build_comparison_summary(metrics_df, cfg.get("comparisons", []))
    comparison_df.to_csv(outdir / "comparison_summary.csv", index=False)

    manifest = {
        **prepared.manifest,
        "n_analyses": int(analysis_registry_df.shape[0]),
        "analysis_blocks": list(analysis_specs_by_block.keys()),
        "binary_endpoints": [endpoint.name for endpoint in binary_endpoints],
        "survival_endpoints": [endpoint.name for endpoint in survival_endpoints],
        "selection_scope": selection_scope,
        "multicollinearity": collinearity_cfg,
        "logistic_selection_rule": logistic_selection_rule,
        "cox_selection_rule": cox_selection_rule,
        "logistic_force_in_note": str(
            logistic_cfg.get(
                "force_in_note",
                "sklearn elastic-net logistic penalizes clinical and immune predictors together in this implementation.",
            )
        ),
        "cox_clinical_penalty_factor": cox_clinical_penalty_factor,
    }
    write_json(manifest, outdir / "run_manifest.json")

    log_step(f"Endpoint modeling completed. Output written to: {outdir}")
    return outdir


def main(config_path: str) -> None:
    run_endpoint_modeling(config_path)


def cli_main() -> None:
    # Keep the CLI intentionally thin; all real work lives in run_endpoint_modeling().
    parser = argparse.ArgumentParser(description="Endpoint modeling for clinical metadata and patient-level immune features.")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    args = parser.parse_args()
    main(args.config)


if __name__ == "__main__":
    cli_main()
