#!/usr/bin/env python

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SRC_ROOT = Path(__file__).resolve().parents[1]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from version1.subtype.clustering import score_patient_clustering, score_patient_clustering_fixed_k
from version1.subtype.io_utils import ensure_outdir, load_run_outputs, read_yaml, write_json


def log_step(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def _resolve_selected_features(outputs: dict) -> list[str]:
    payload = outputs.get("selected_features")
    if isinstance(payload, dict) and isinstance(payload.get("selected_features"), list):
        return [str(x) for x in payload["selected_features"]]

    metrics = outputs.get("metrics")
    if isinstance(metrics, dict) and isinstance(metrics.get("selected_feature_names"), list):
        return [str(x) for x in metrics["selected_feature_names"]]

    raise ValueError("Could not find selected features in selected_features.json or metrics.json")


def _resolve_config(outputs: dict, config_path: str | None) -> dict:
    if config_path:
        return read_yaml(config_path)
    if outputs.get("config") is not None:
        return outputs["config"]
    raise ValueError("Could not find config.yaml in the run directory; please pass --config explicitly.")


def run_permutation_importance(
    run_dir: Path,
    n_permutations: int,
    seed: int,
    config_path: str | None = None,
    outdir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    outputs = load_run_outputs(run_dir)
    patient_features = outputs.get("patient_features")
    if patient_features is None:
        raise FileNotFoundError(f"patient_features.csv was not found under {run_dir}")

    cfg = _resolve_config(outputs, config_path)
    selected_features = _resolve_selected_features(outputs)
    if len(selected_features) == 0:
        raise ValueError("No selected features were found for clustering.")

    missing = [f for f in selected_features if f not in patient_features.columns]
    if missing:
        raise ValueError(f"Selected features are missing from patient_features.csv: {missing[:10]}")

    outdir = ensure_outdir(outdir or (run_dir / "permutation_feature_importance"))

    log_step("Computing baseline clustering score")
    base_cluster_df, base_silhouette_df, base_metrics = score_patient_clustering(
        patient_features,
        selected_features,
        cfg,
    )
    base_silhouette_free_k = base_metrics.get("best_silhouette", np.nan)
    base_best_k = int(base_metrics.get("best_k", 1))
    if pd.isna(base_silhouette_free_k):
        raise ValueError(
            "Baseline best_silhouette is NaN. Permutation importance requires a valid silhouette score."
        )
    _, base_silhouette_fixed_df, base_fixed_metrics = score_patient_clustering_fixed_k(
        patient_features,
        selected_features,
        cfg,
        fixed_k=base_best_k,
    )
    base_silhouette_fixed_k = base_fixed_metrics.get("best_silhouette", np.nan)

    rng = np.random.default_rng(seed)
    rows = []

    for idx, feature in enumerate(selected_features, start=1):
        log_step(f"Permuting feature {idx}/{len(selected_features)}: {feature}")
        original_values = patient_features[feature].to_numpy(copy=True)

        for perm_idx in range(1, n_permutations + 1):
            permuted_df = patient_features.copy()
            permuted_df[feature] = rng.permutation(original_values)

            _, _, perm_metrics_free = score_patient_clustering(
                permuted_df,
                selected_features,
                cfg,
            )
            _, _, perm_metrics_fixed = score_patient_clustering_fixed_k(
                permuted_df,
                selected_features,
                cfg,
                fixed_k=base_best_k,
            )

            rows.append({
                "feature": feature,
                "permutation": perm_idx,
                "base_best_k": int(base_best_k),
                "base_silhouette_free_k": float(base_silhouette_free_k),
                "base_silhouette_fixed_k": float(base_silhouette_fixed_k) if pd.notna(base_silhouette_fixed_k) else np.nan,
                "perm_silhouette_free_k": float(perm_metrics_free.get("best_silhouette", np.nan)),
                "perm_best_k_free": int(perm_metrics_free.get("best_k", 1)),
                "perm_silhouette_fixed_k": float(perm_metrics_fixed.get("best_silhouette", np.nan)),
            })

    long_df = pd.DataFrame(rows)
    if long_df.empty:
        raise ValueError("No permutation results were generated.")

    summary_df = (
        long_df.groupby("feature", as_index=False)
        .agg(
            n_permutations=("permutation", "size"),
            base_best_k=("base_best_k", "first"),
            base_silhouette_free_k=("base_silhouette_free_k", "first"),
            base_silhouette_fixed_k=("base_silhouette_fixed_k", "first"),
            mean_perm_silhouette_free_k=("perm_silhouette_free_k", "mean"),
            std_perm_silhouette_free_k=("perm_silhouette_free_k", "std"),
            min_perm_silhouette_free_k=("perm_silhouette_free_k", "min"),
            max_perm_silhouette_free_k=("perm_silhouette_free_k", "max"),
            mean_perm_best_k_free=("perm_best_k_free", "mean"),
            mean_perm_silhouette_fixed_k=("perm_silhouette_fixed_k", "mean"),
            std_perm_silhouette_fixed_k=("perm_silhouette_fixed_k", "std"),
            min_perm_silhouette_fixed_k=("perm_silhouette_fixed_k", "min"),
            max_perm_silhouette_fixed_k=("perm_silhouette_fixed_k", "max"),
        )
    )
    summary_df["importance_free_k"] = summary_df["base_silhouette_free_k"] - summary_df["mean_perm_silhouette_free_k"]
    summary_df["importance_fixed_k"] = summary_df["base_silhouette_fixed_k"] - summary_df["mean_perm_silhouette_fixed_k"]
    summary_df["abs_importance_free_k"] = summary_df["importance_free_k"].abs()
    summary_df["abs_importance_fixed_k"] = summary_df["importance_fixed_k"].abs()
    summary_df = summary_df.sort_values("importance_fixed_k", ascending=False).reset_index(drop=True)
    summary_df.insert(0, "rank", np.arange(1, summary_df.shape[0] + 1))

    long_df.to_csv(outdir / "permutation_importance_long.csv", index=False)
    summary_df.to_csv(outdir / "permutation_importance_summary.csv", index=False)
    base_cluster_df.to_csv(outdir / "baseline_clusters.csv", index=False)
    base_silhouette_df.to_csv(outdir / "baseline_silhouette_by_k.csv", index=False)
    base_silhouette_fixed_df.to_csv(outdir / "baseline_silhouette_fixed_k.csv", index=False)
    write_json(
        {
            "run_dir": str(run_dir),
            "n_permutations": int(n_permutations),
            "seed": int(seed),
            "base_best_k": int(base_best_k),
            "base_best_silhouette_free_k": float(base_silhouette_free_k),
            "base_best_silhouette_fixed_k": float(base_silhouette_fixed_k) if pd.notna(base_silhouette_fixed_k) else None,
            "selected_features": selected_features,
        },
        outdir / "permutation_importance_meta.json",
    )

    return summary_df, long_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Permutation-based feature importance for patient clustering.")
    parser.add_argument("--run-dir", required=True, help="Subtype run directory containing patient_features.csv")
    parser.add_argument("--config", help="Optional config path. Defaults to <run-dir>/config.yaml")
    parser.add_argument("--n-permutations", type=int, default=20, help="Number of permutations per feature")
    parser.add_argument("--seed", type=int, default=123, help="Random seed for feature permutation")
    parser.add_argument(
        "--outdir",
        help="Optional output directory. Defaults to <run-dir>/permutation_feature_importance",
    )
    args = parser.parse_args()

    summary_df, _ = run_permutation_importance(
        run_dir=Path(args.run_dir),
        n_permutations=int(args.n_permutations),
        seed=int(args.seed),
        config_path=args.config,
        outdir=Path(args.outdir) if args.outdir else None,
    )

    print(summary_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
