#!/usr/bin/env python

from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact

SRC_ROOT = Path(__file__).resolve().parents[1]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from version1.subtype.io_utils import ensure_outdir, load_run_outputs, read_yaml, write_json


def log_step(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def require_optional_dependencies():
    try:
        import matplotlib.pyplot as plt
        from lifelines import CoxPHFitter, KaplanMeierFitter
        from lifelines.plotting import add_at_risk_counts
        from lifelines.statistics import multivariate_logrank_test
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "cluster_clinical_association.py requires optional packages "
            "`openpyxl`, `lifelines`, and `matplotlib`. "
            "Install them in the active environment first."
        ) from e
    return plt, CoxPHFitter, KaplanMeierFitter, add_at_risk_counts, multivariate_logrank_test


def resolve_config(outputs: dict, config_path: str | None) -> dict:
    if config_path:
        return read_yaml(config_path)
    if outputs.get("config") is not None:
        return outputs["config"]
    raise ValueError("Could not find config.yaml in the run directory; please pass --config explicitly.")


def read_clinical_excel(path: str | Path, sheet_name) -> pd.DataFrame:
    try:
        return pd.read_excel(path, sheet_name=sheet_name, engine="openpyxl")
    except ImportError as e:
        raise ImportError(
            "Reading the clinical Excel file requires `openpyxl`. "
            "Install it in the active environment first."
        ) from e


def as_flag(series: pd.Series, positive_values: list) -> pd.Series:
    positive = {str(x) for x in positive_values}
    return series.astype(str).isin(positive).astype(int)


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    try:
        from statsmodels.stats.multitest import multipletests
    except ModuleNotFoundError:
        return pd.Series(np.nan, index=pvalues.index)

    valid = pvalues.notna()
    out = pd.Series(np.nan, index=pvalues.index, dtype=float)
    if valid.any():
        out.loc[valid] = multipletests(pvalues.loc[valid], method="fdr_bh")[1]
    return out


def prepare_merged_data(run_dir: Path, cfg: dict) -> tuple[pd.DataFrame, dict]:
    outputs = load_run_outputs(run_dir)
    patient_features = outputs.get("patient_features")
    cluster_df = outputs.get("patient_clusters")
    if patient_features is None or cluster_df is None:
        raise FileNotFoundError("Run directory must contain patient_features.csv and patient_clusters.csv")

    ca_cfg = cfg.get("clinical_association", {})
    excel_path = ca_cfg.get("excel_path")
    if not excel_path:
        raise ValueError("clinical_association.excel_path is required in the config.")

    run_patient_key = ca_cfg.get("run_patient_key", "patient_id")
    clinical_join_key = ca_cfg.get("clinical_join_key", "sample_id")

    clinical_df = read_clinical_excel(excel_path, ca_cfg.get("sheet_name", 0))
    if clinical_join_key not in clinical_df.columns:
        raise ValueError(f"Clinical Excel is missing join column: {clinical_join_key}")
    if run_patient_key not in patient_features.columns:
        raise ValueError(f"patient_features.csv is missing join column: {run_patient_key}")

    run_meta = patient_features[["analysis_id", run_patient_key]].copy()
    if "source" in patient_features.columns:
        run_meta["source"] = patient_features["source"].values
    if "batch" in patient_features.columns:
        run_meta["batch"] = patient_features["batch"].values

    run_meta[run_patient_key] = run_meta[run_patient_key].astype(str)
    clinical_df = clinical_df.copy()
    clinical_df[clinical_join_key] = clinical_df[clinical_join_key].astype(str)

    merged = (
        cluster_df.merge(run_meta, on="analysis_id", how="left")
        .merge(clinical_df, left_on=run_patient_key, right_on=clinical_join_key, how="inner")
    )

    meta = {
        "n_clustered_patients": int(cluster_df["analysis_id"].nunique()),
        "n_patients_with_clinical_match": int(merged["analysis_id"].nunique()),
        "clinical_excel_path": str(excel_path),
        "clinical_join_key": clinical_join_key,
        "run_patient_key": run_patient_key,
    }
    return merged, meta


def run_fisher_tests(dat: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    ca_cfg = cfg["clinical_association"]
    cluster_col = ca_cfg.get("cluster_column", "patient_cluster")
    benefit_col = ca_cfg.get("benefit_column", "PD_Event")
    benefit_positive_values = ca_cfg.get("benefit_positive_values", [0])

    fisher_df = dat[[cluster_col, benefit_col]].copy()
    fisher_df["benefit_flag"] = as_flag(fisher_df[benefit_col], benefit_positive_values)
    fisher_df = fisher_df.dropna(subset=[cluster_col])

    rows = []
    for cluster in sorted(fisher_df[cluster_col].dropna().unique()):
        in_cluster = fisher_df[cluster_col] == cluster
        benefit = fisher_df["benefit_flag"] == 1

        a = int((in_cluster & benefit).sum())
        b = int((in_cluster & ~benefit).sum())
        c = int((~in_cluster & benefit).sum())
        d = int((~in_cluster & ~benefit).sum())

        odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        rows.append({
            "cluster": int(cluster),
            "n_in_cluster": int(in_cluster.sum()),
            "n_benefit_in_cluster": a,
            "n_nonbenefit_in_cluster": b,
            "n_benefit_out_cluster": c,
            "n_nonbenefit_out_cluster": d,
            "odds_ratio": float(odds_ratio),
            "pvalue": float(pvalue),
        })

    out = pd.DataFrame(rows).sort_values("pvalue", ascending=True).reset_index(drop=True)
    if not out.empty:
        out["pvalue_fdr_bh"] = bh_adjust(out["pvalue"])
    return out


def run_kaplan_meier(dat: pd.DataFrame, cfg: dict, outdir: Path):
    plt, _, KaplanMeierFitter, add_at_risk_counts, multivariate_logrank_test = require_optional_dependencies()

    ca_cfg = cfg["clinical_association"]
    cluster_col = ca_cfg.get("cluster_column", "patient_cluster")
    time_col = ca_cfg.get("pfs_time_column", "PFS (Days)")
    event_col = ca_cfg.get("pfs_event_column", "PD_Event")
    event_positive_values = ca_cfg.get("pfs_event_positive_values", [1])

    km_df = dat[[cluster_col, time_col, event_col]].copy()
    km_df["PFS_time"] = pd.to_numeric(km_df[time_col], errors="coerce")
    km_df["PFS_event"] = as_flag(km_df[event_col], event_positive_values)
    km_df = km_df.dropna(subset=[cluster_col, "PFS_time"])

    if km_df.empty or km_df[cluster_col].nunique() < 2:
        return pd.DataFrame(), {"logrank_pvalue": np.nan, "n_patients": int(km_df.shape[0])}

    result = multivariate_logrank_test(
        event_durations=km_df["PFS_time"],
        groups=km_df[cluster_col],
        event_observed=km_df["PFS_event"],
    )

    fig, ax = plt.subplots(figsize=(7, 6))
    km_fitters = []
    rows = []
    for cluster in sorted(km_df[cluster_col].unique()):
        sub = km_df[km_df[cluster_col] == cluster]
        kmf = KaplanMeierFitter(label=f"C{int(cluster)}")
        kmf.fit(sub["PFS_time"], event_observed=sub["PFS_event"])
        kmf.plot_survival_function(ax=ax, ci_show=False)
        km_fitters.append(kmf)
        rows.append({
            "cluster": int(cluster),
            "n_patients": int(sub.shape[0]),
            "n_events": int(sub["PFS_event"].sum()),
            "median_pfs_time": float(kmf.median_survival_time_),
        })

    add_at_risk_counts(*km_fitters, ax=ax)
    ax.set_xlabel("PFS time (days)")
    ax.set_ylabel("Survival probability")
    ax.set_title(f"PFS by cluster (log-rank p={result.p_value:.3g})")
    fig.tight_layout()
    fig.savefig(outdir / "kaplan_meier_pfs_by_cluster.png", dpi=200)
    plt.close(fig)

    return pd.DataFrame(rows), {
        "logrank_test_statistic": float(result.test_statistic),
        "logrank_pvalue": float(result.p_value),
        "n_patients": int(km_df.shape[0]),
    }


def build_cox_design_matrix(dat: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, list[str]]:
    ca_cfg = cfg["clinical_association"]
    cluster_col = ca_cfg.get("cluster_column", "patient_cluster")
    time_col = ca_cfg.get("pfs_time_column", "PFS (Days)")
    event_col = ca_cfg.get("pfs_event_column", "PD_Event")
    event_positive_values = ca_cfg.get("pfs_event_positive_values", [1])
    covariates = [str(x) for x in ca_cfg.get("cox_covariates", [])]
    categorical_covariates = {str(x) for x in ca_cfg.get("categorical_covariates", [])}

    design_parts = []
    core = pd.DataFrame({
        "PFS_time": pd.to_numeric(dat[time_col], errors="coerce"),
        "PFS_event": as_flag(dat[event_col], event_positive_values),
    })
    design_parts.append(core)

    cluster_dummies = pd.get_dummies(
        dat[cluster_col].astype("Int64").astype(str),
        prefix="cluster",
        drop_first=True,
    ).astype(int)
    cluster_cols = cluster_dummies.columns.tolist()
    design_parts.append(cluster_dummies)

    for cov in covariates:
        if cov not in dat.columns:
            continue
        if cov in categorical_covariates:
            dummies = pd.get_dummies(
                dat[cov].astype("string").fillna("NA"),
                prefix=sanitize_name(cov),
                drop_first=True,
            ).astype(int)
            if dummies.shape[1] > 0:
                design_parts.append(dummies)
        else:
            design_parts.append(
                pd.DataFrame({sanitize_name(cov): pd.to_numeric(dat[cov], errors="coerce")})
            )

    cox_df = pd.concat(design_parts, axis=1)
    cox_df = cox_df.replace([np.inf, -np.inf], np.nan).dropna()
    return cox_df, cluster_cols


def run_adjusted_cox(dat: pd.DataFrame, cfg: dict):
    _, CoxPHFitter, _, _, _ = require_optional_dependencies()

    cox_df, cluster_cols = build_cox_design_matrix(dat, cfg)
    if cox_df.empty or len(cluster_cols) == 0:
        empty = pd.DataFrame()
        return empty, empty, {"cluster_global_lr_pvalue": np.nan, "n_patients": int(cox_df.shape[0])}

    full = CoxPHFitter()
    full.fit(cox_df, duration_col="PFS_time", event_col="PFS_event")
    full_summary = full.summary.reset_index().rename(columns={"covariate": "term"})

    reduced_df = cox_df.drop(columns=cluster_cols)
    reduced = CoxPHFitter()
    reduced.fit(reduced_df, duration_col="PFS_time", event_col="PFS_event")

    lr_stat = 2.0 * (full.log_likelihood_ - reduced.log_likelihood_)
    lr_df = int(len(cluster_cols))
    lr_pvalue = float(chi2.sf(lr_stat, lr_df))

    cluster_summary = full_summary[full_summary["term"].isin(cluster_cols)].copy()
    cluster_summary["hazard_ratio"] = np.exp(cluster_summary["coef"])

    meta = {
        "n_patients": int(cox_df.shape[0]),
        "n_cluster_terms": int(len(cluster_cols)),
        "cluster_global_lr_statistic": float(lr_stat),
        "cluster_global_lr_df": lr_df,
        "cluster_global_lr_pvalue": lr_pvalue,
    }
    return cluster_summary, full_summary, meta


def main() -> None:
    parser = argparse.ArgumentParser(description="Cluster-clinical association analysis for subtype runs.")
    parser.add_argument("--run-dir", required=True, help="Subtype run directory")
    parser.add_argument("--config", help="Optional config path. Defaults to <run-dir>/config.yaml")
    parser.add_argument(
        "--outdir",
        help="Optional output directory. Defaults to <run-dir>/clinical_association",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    outputs = load_run_outputs(run_dir)
    cfg = resolve_config(outputs, args.config)
    outdir = ensure_outdir(Path(args.outdir) if args.outdir else run_dir / "clinical_association")

    log_step("Loading cluster outputs and clinical Excel")
    merged, merge_meta = prepare_merged_data(run_dir, cfg)
    merged.to_csv(outdir / "cluster_clinical_merged.csv", index=False)

    log_step("Running Fisher enrichment tests for clinical benefit")
    fisher_df = run_fisher_tests(merged, cfg)
    fisher_df.to_csv(outdir / "cluster_benefit_fisher.csv", index=False)

    log_step("Running Kaplan-Meier analysis for PFS")
    km_df, km_meta = run_kaplan_meier(merged, cfg, outdir)
    km_df.to_csv(outdir / "kaplan_meier_cluster_summary.csv", index=False)

    log_step("Running adjusted Cox model")
    cox_cluster_df, cox_full_df, cox_meta = run_adjusted_cox(merged, cfg)
    cox_cluster_df.to_csv(outdir / "adjusted_cox_cluster_terms.csv", index=False)
    cox_full_df.to_csv(outdir / "adjusted_cox_full_summary.csv", index=False)

    write_json(
        {
            "merge": merge_meta,
            "kaplan_meier": km_meta,
            "cox": cox_meta,
        },
        outdir / "clinical_association_meta.json",
    )

    print(fisher_df.to_string(index=False) if not fisher_df.empty else "No Fisher results.")
    if not cox_cluster_df.empty:
        print()
        print(cox_cluster_df.to_string(index=False))


if __name__ == "__main__":
    main()
