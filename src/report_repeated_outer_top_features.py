#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SUMMARY_BASENAME = "repeated_top_features"


def _outer_validation_dir(run_dir: Path) -> Path:
    outer_dir = run_dir / "outer_search_validation"
    if not outer_dir.exists():
        raise FileNotFoundError(f"Missing outer_search_validation directory: {outer_dir}")
    return outer_dir


def _discover_endpoint_slugs(outer_dir: Path) -> list[str]:
    return sorted([p.name for p in outer_dir.iterdir() if p.is_dir()])


def _feature_catalog(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "feature_catalog_resolved.csv"
    if not path.exists():
        return pd.DataFrame(columns=["feature_name", "family", "family_cap_group", "feature_level", "description"])
    keep_cols = ["feature_name", "family", "family_cap_group", "feature_level", "description"]
    df = pd.read_csv(path)
    cols = [col for col in keep_cols if col in df.columns]
    return df[cols].copy()


def _group_key_from_columns(df: pd.DataFrame, feature_col: str) -> pd.Series:
    if "family_cap_group" in df.columns:
        family_cap_group = df["family_cap_group"]
        return family_cap_group.where(family_cap_group.notna() & (family_cap_group.astype(str) != ""), df[feature_col])
    return df[feature_col]


def _load_outer_selected(outer_dir: Path) -> pd.DataFrame:
    path = outer_dir / "outer_selected_candidates.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _load_outer_fold_metrics(outer_dir: Path) -> pd.DataFrame:
    path = outer_dir / "outer_fold_metrics.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _fold_stage2_paths(outer_dir: Path, endpoint_slug: str) -> list[Path]:
    endpoint_dir = outer_dir / endpoint_slug
    if not endpoint_dir.exists():
        raise FileNotFoundError(f"Missing endpoint directory: {endpoint_dir}")
    return sorted(endpoint_dir.glob("fold_*/stage2_single_feature_scan.csv"))


def _build_outer_eval_lookup(
    outer_selected: pd.DataFrame, outer_fold_metrics: pd.DataFrame, feature_catalog: pd.DataFrame
) -> pd.DataFrame:
    if outer_selected.empty or outer_fold_metrics.empty:
        return pd.DataFrame(
            columns=[
                "endpoint_slug",
                "fold",
                "selected_group_key",
                "outer_selected_feature_name",
                "outer_selected_roc_auc",
                "outer_selected_auprc",
                "outer_selected_brier",
                "outer_selected_cindex",
                "outer_delta_roc_auc",
                "outer_delta_auprc",
                "outer_delta_brier",
                "outer_delta_cindex",
            ]
        )

    selected_rows = outer_fold_metrics[outer_fold_metrics["model_role"] == "selected"].copy()
    baseline_rows = outer_fold_metrics[outer_fold_metrics["model_role"] == "baseline"].copy()
    merge_cols = ["endpoint_slug", "fold"]
    baseline_metric_cols = [col for col in ["roc_auc", "auprc", "brier", "cindex"] if col in baseline_rows.columns]
    selected_metric_cols = [col for col in ["roc_auc", "auprc", "brier", "cindex"] if col in selected_rows.columns]

    selected_rows = selected_rows[merge_cols + selected_metric_cols].rename(
        columns={col: f"outer_selected_{col}" for col in selected_metric_cols}
    )
    baseline_rows = baseline_rows[merge_cols + baseline_metric_cols].rename(
        columns={col: f"outer_baseline_{col}" for col in baseline_metric_cols}
    )
    merged = outer_selected.merge(selected_rows, on=merge_cols, how="left").merge(baseline_rows, on=merge_cols, how="left")
    if "feature_names" in merged.columns:
        merged = merged.rename(columns={"feature_names": "outer_selected_feature_name"})

    if "outer_selected_feature_name" in merged.columns and not feature_catalog.empty:
        group_lookup = feature_catalog[["feature_name", "family_cap_group"]].drop_duplicates().rename(
            columns={"feature_name": "outer_selected_feature_name"}
        )
        merged = merged.merge(group_lookup, on="outer_selected_feature_name", how="left")
        merged["selected_group_key"] = _group_key_from_columns(merged, "outer_selected_feature_name")
        merged = merged.drop(columns=[col for col in ["family_cap_group"] if col in merged.columns])
    elif "outer_selected_feature_name" in merged.columns:
        merged["selected_group_key"] = merged["outer_selected_feature_name"]
    else:
        merged["selected_group_key"] = np.nan

    for metric in ["roc_auc", "auprc", "brier", "cindex"]:
        sel_col = f"outer_selected_{metric}"
        base_col = f"outer_baseline_{metric}"
        delta_col = f"outer_delta_{metric}"
        if sel_col in merged.columns and base_col in merged.columns:
            merged[delta_col] = merged[sel_col] - merged[base_col]
        else:
            merged[delta_col] = np.nan
    keep_cols = [
        "endpoint_slug",
        "fold",
        "selected_group_key",
        "outer_selected_feature_name",
        "outer_selected_roc_auc",
        "outer_selected_auprc",
        "outer_selected_brier",
        "outer_selected_cindex",
        "outer_delta_roc_auc",
        "outer_delta_auprc",
        "outer_delta_brier",
        "outer_delta_cindex",
    ]
    return merged[[col for col in keep_cols if col in merged.columns]].copy()


def _metric_columns(df: pd.DataFrame) -> list[str]:
    preferred = [
        "ranking_metric",
        "roc_auc_mean",
        "delta_roc_auc_mean",
        "auprc_mean",
        "delta_auprc_mean",
        "brier_mean",
        "delta_brier_mean",
        "cindex_mean",
        "delta_cindex_mean",
        "full_inner_score",
    ]
    return [col for col in preferred if col in df.columns]


def _fold_label(fold: int) -> str:
    return f"{int(fold):02d}"


def _representative_features(details_df: pd.DataFrame) -> pd.DataFrame:
    rep_df = (
        details_df.groupby(["endpoint_slug", "group_key", "feature_name"], as_index=False)
        .agg(
            n_folds_for_feature=("fold", "nunique"),
            mean_group_rank_for_feature=("group_rank_in_fold", "mean"),
            mean_ranking_metric_for_feature=("ranking_metric", "mean"),
            family=("family", "first"),
            family_cap_group=("family_cap_group", "first"),
            feature_level=("feature_level", "first"),
            description=("description", "first"),
        )
        .sort_values(
            ["endpoint_slug", "group_key", "mean_group_rank_for_feature", "n_folds_for_feature", "mean_ranking_metric_for_feature"],
            ascending=[True, True, True, False, False],
        )
        .drop_duplicates(["endpoint_slug", "group_key"], keep="first")
        .rename(columns={"feature_name": "representative_feature_name"})
    )
    features_seen_df = (
        details_df.groupby(["endpoint_slug", "group_key"])["feature_name"]
        .agg(lambda s: "|".join(sorted(pd.unique(s))))
        .reset_index(name="features_seen_in_group")
    )
    return rep_df.merge(features_seen_df, on=["endpoint_slug", "group_key"], how="left")


def _pivot_by_fold(details_df: pd.DataFrame, value_col: str, prefix: str) -> pd.DataFrame:
    if value_col not in details_df.columns:
        return pd.DataFrame(columns=["endpoint_slug", "group_key"])
    wide = details_df.pivot_table(index=["endpoint_slug", "group_key"], columns="fold", values=value_col, aggfunc="first")
    if wide.empty:
        return pd.DataFrame(columns=["endpoint_slug", "group_key"])
    wide = wide.rename(columns={fold: f"{prefix}_fold_{_fold_label(fold)}" for fold in wide.columns}).reset_index()
    return wide


def collect_top_features(run_dir: Path, endpoint_slugs: list[str], top_n: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    outer_dir = _outer_validation_dir(run_dir)
    feature_catalog = _feature_catalog(run_dir)
    outer_selected = _load_outer_selected(outer_dir)
    outer_fold_metrics = _load_outer_fold_metrics(outer_dir)
    outer_eval_lookup = _build_outer_eval_lookup(outer_selected, outer_fold_metrics, feature_catalog)

    detail_rows: list[dict[str, object]] = []

    for endpoint_slug in endpoint_slugs:
        for stage2_path in _fold_stage2_paths(outer_dir, endpoint_slug):
            fold_name = stage2_path.parent.name
            fold_num = int(fold_name.split("_")[-1])
            df = pd.read_csv(stage2_path)
            if "status" in df.columns:
                df = df[df["status"] == "ok"].copy()
            if df.empty:
                continue
            if "endpoint_slug" in df.columns:
                df = df[df["endpoint_slug"] == endpoint_slug].copy()
            if "ranking_metric" not in df.columns:
                raise ValueError(f"Missing ranking_metric in {stage2_path}")
            df = df.sort_values(["ranking_metric", "full_inner_score"], ascending=[False, False]).copy()
            df["raw_rank_in_fold"] = range(1, len(df) + 1)
            df["endpoint_slug"] = endpoint_slug
            df["fold"] = fold_num
            if not feature_catalog.empty:
                df = df.merge(feature_catalog, on="feature_name", how="left")
            df["group_key"] = _group_key_from_columns(df, "feature_name")
            # Collapse highly related score variants to one representative per family_cap_group within each fold.
            df = df.drop_duplicates(subset=["group_key"], keep="first").copy()
            df["group_rank_in_fold"] = range(1, len(df) + 1)
            df["in_top_n"] = df["group_rank_in_fold"] <= top_n
            metric_cols = _metric_columns(df)
            keep_cols = [
                "endpoint_slug",
                "fold",
                "feature_name",
                "group_key",
                "family",
                "family_cap_group",
                "feature_level",
                "description",
                "raw_rank_in_fold",
                "group_rank_in_fold",
                "in_top_n",
            ] + metric_cols
            fold_df = df[[col for col in keep_cols if col in df.columns]].copy()
            detail_rows.extend(fold_df.to_dict(orient="records"))

    if not detail_rows:
        empty = pd.DataFrame()
        return empty, empty

    details_df = pd.DataFrame(detail_rows)

    if not outer_eval_lookup.empty:
        details_df = details_df.merge(
            outer_eval_lookup,
            left_on=["endpoint_slug", "fold", "group_key"],
            right_on=["endpoint_slug", "fold", "selected_group_key"],
            how="left",
        )
        selected_metric_cols = [
            col
            for col in [
                "outer_selected_roc_auc",
                "outer_selected_auprc",
                "outer_selected_brier",
                "outer_selected_cindex",
            ]
            if col in details_df.columns
        ]
        if selected_metric_cols:
            details_df["selected_in_outer_fold"] = details_df[selected_metric_cols].notna().any(axis=1)
        else:
            details_df["selected_in_outer_fold"] = False
        details_df = details_df.drop(
            columns=[col for col in ["selected_group_key"] if col in details_df.columns]
        )
    else:
        details_df["selected_in_outer_fold"] = False

    representative_df = _representative_features(details_df)

    agg_spec: dict[str, tuple[str, str]] = {
        "n_folds_total": ("fold", "nunique"),
        "folds_seen": ("fold", lambda s: ",".join(map(str, sorted(set(int(x) for x in s))))),
        "mean_raw_rank": ("raw_rank_in_fold", "mean"),
        "mean_group_rank": ("group_rank_in_fold", "mean"),
        "median_group_rank": ("group_rank_in_fold", "median"),
        "best_group_rank": ("group_rank_in_fold", "min"),
    }
    for col in _metric_columns(details_df):
        agg_spec[f"mean_{col}"] = (col, "mean")
    if "selected_in_outer_fold" in details_df.columns:
        agg_spec["selected_in_outer_n_folds"] = ("selected_in_outer_fold", "sum")
    for col in [
        "outer_selected_roc_auc",
        "outer_selected_auprc",
        "outer_selected_brier",
        "outer_selected_cindex",
        "outer_delta_roc_auc",
        "outer_delta_auprc",
        "outer_delta_brier",
        "outer_delta_cindex",
    ]:
        if col in details_df.columns:
            agg_spec[f"mean_{col}"] = (col, "mean")

    summary_df = (
        details_df.groupby(["endpoint_slug", "group_key"], as_index=False)
        .agg(**agg_spec)
        .sort_values(["endpoint_slug", "mean_group_rank", "mean_ranking_metric"], ascending=[True, True, False])
        .reset_index(drop=True)
    )

    top_n_hits = (
        details_df[details_df["in_top_n"]]
        .groupby(["endpoint_slug", "group_key"], as_index=False)
        .agg(
            n_folds_in_top_n=("fold", "nunique"),
            folds_in_top_n=("fold", lambda s: ",".join(map(str, sorted(set(int(x) for x in s))))),
        )
    )
    summary_df = summary_df.merge(top_n_hits, on=["endpoint_slug", "group_key"], how="left")
    summary_df["n_folds_in_top_n"] = summary_df["n_folds_in_top_n"].fillna(0).astype(int)
    summary_df["folds_in_top_n"] = summary_df["folds_in_top_n"].fillna("")
    summary_df["top_n_reference"] = top_n

    total_folds_df = (
        details_df.groupby("endpoint_slug", as_index=False)["fold"]
        .nunique()
        .rename(columns={"fold": "total_outer_folds"})
    )
    summary_df = summary_df.merge(total_folds_df, on="endpoint_slug", how="left")
    summary_df = summary_df.merge(representative_df, on=["endpoint_slug", "group_key"], how="left")

    fold_level_wide = [
        _pivot_by_fold(details_df, "raw_rank_in_fold", "raw_rank"),
        _pivot_by_fold(details_df, "group_rank_in_fold", "group_rank"),
        _pivot_by_fold(details_df, "ranking_metric", "ranking_metric"),
        _pivot_by_fold(details_df, "roc_auc_mean", "roc_auc_mean"),
        _pivot_by_fold(details_df, "delta_roc_auc_mean", "delta_roc_auc_mean"),
        _pivot_by_fold(details_df, "outer_selected_roc_auc", "outer_selected_roc_auc"),
        _pivot_by_fold(details_df, "outer_delta_roc_auc", "outer_delta_roc_auc"),
    ]
    for wide_df in fold_level_wide:
        if len(wide_df.columns) > 2:
            summary_df = summary_df.merge(wide_df, on=["endpoint_slug", "group_key"], how="left")

    ordered_cols = [
        "endpoint_slug",
        "group_key",
        "representative_feature_name",
        "features_seen_in_group",
        "family",
        "family_cap_group",
        "feature_level",
        "description",
        "top_n_reference",
        "n_folds_in_top_n",
        "total_outer_folds",
        "folds_in_top_n",
        "n_folds_total",
        "folds_seen",
        "mean_raw_rank",
        "mean_group_rank",
        "median_group_rank",
        "best_group_rank",
        "mean_ranking_metric",
        "mean_roc_auc_mean",
        "mean_delta_roc_auc_mean",
        "mean_auprc_mean",
        "mean_delta_auprc_mean",
        "mean_brier_mean",
        "mean_delta_brier_mean",
        "mean_cindex_mean",
        "mean_delta_cindex_mean",
        "selected_in_outer_n_folds",
        "mean_outer_selected_roc_auc",
        "mean_outer_delta_roc_auc",
        "mean_outer_selected_auprc",
        "mean_outer_delta_auprc",
        "mean_outer_selected_brier",
        "mean_outer_delta_brier",
        "mean_outer_selected_cindex",
        "mean_outer_delta_cindex",
    ]
    ordered_cols.extend(
        sorted(
            [
                col
                for col in summary_df.columns
                if col.startswith(
                    (
                        "raw_rank_fold_",
                        "group_rank_fold_",
                        "ranking_metric_fold_",
                        "roc_auc_mean_fold_",
                        "delta_roc_auc_mean_fold_",
                        "outer_selected_roc_auc_fold_",
                        "outer_delta_roc_auc_fold_",
                    )
                )
            ]
        )
    )
    summary_df = summary_df[[col for col in ordered_cols if col in summary_df.columns]]

    detail_order = [
        "endpoint_slug",
        "fold",
        "group_key",
        "feature_name",
        "outer_selected_feature_name",
        "family",
        "family_cap_group",
        "feature_level",
        "description",
        "raw_rank_in_fold",
        "group_rank_in_fold",
        "in_top_n",
        "ranking_metric",
        "roc_auc_mean",
        "delta_roc_auc_mean",
        "auprc_mean",
        "delta_auprc_mean",
        "brier_mean",
        "delta_brier_mean",
        "cindex_mean",
        "delta_cindex_mean",
        "full_inner_score",
        "selected_in_outer_fold",
        "outer_selected_roc_auc",
        "outer_delta_roc_auc",
        "outer_selected_auprc",
        "outer_delta_auprc",
        "outer_selected_brier",
        "outer_delta_brier",
        "outer_selected_cindex",
        "outer_delta_cindex",
    ]
    details_df = details_df[[col for col in detail_order if col in details_df.columns]].sort_values(
        ["endpoint_slug", "fold", "group_rank_in_fold", "feature_name"]
    )
    return summary_df, details_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect features that repeatedly appear in the top-N stage2 outer-validation rankings."
    )
    parser.add_argument("--run-dir", required=True, help="Run directory, e.g. data/.../feature_search_base_v2_single_feature_outer")
    parser.add_argument("--top-n", type=int, default=10, help="Top-N features to collect from each outer fold stage2 scan")
    parser.add_argument(
        "--endpoint-slug",
        action="append",
        dest="endpoint_slugs",
        help="Endpoint slug(s) under outer_search_validation. If omitted, all endpoint directories are used.",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    outer_dir = _outer_validation_dir(run_dir)
    endpoint_slugs = args.endpoint_slugs or _discover_endpoint_slugs(outer_dir)
    if not endpoint_slugs:
        raise ValueError(f"No endpoint directories found under {outer_dir}")
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")

    summary_df, details_df = collect_top_features(run_dir, endpoint_slugs, args.top_n)

    summary_path = outer_dir / f"{SUMMARY_BASENAME}__top{args.top_n}.csv"
    detail_path = outer_dir / f"{SUMMARY_BASENAME}__top{args.top_n}__details.csv"
    summary_df.to_csv(summary_path, index=False)
    details_df.to_csv(detail_path, index=False)

    print(f"Wrote summary: {summary_path}")
    print(f"Wrote details: {detail_path}")
    if summary_df.empty:
        print("No top-N features found.")
    else:
        print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
