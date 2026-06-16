from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Iterable, Optional, Dict, Any
import logging

import scanpy as sc


def annotate_gene_sets(adata, mt_prefix: str = "MT-", ribo_prefixes=("RPS", "RPL"), hb_prefix=("HBA", "HBB")):
    """
    Adds boolean columns in adata.var for mitochondrial/ribosomal/hemoglobin genes.
    """
    genes = adata.var_names.astype(str)
    adata.var["mt"] = genes.str.upper().str.startswith(mt_prefix.upper())
    adata.var["ribo"] = genes.str.upper().str.startswith(tuple(p.upper() for p in ribo_prefixes))
    adata.var["hb"] = genes.str.upper().str.startswith(tuple(p.upper() for p in hb_prefix))
    return adata


def compute_qc_metrics(adata, logger: Optional[logging.Logger] = None):
    """
    Computes per-cell QC metrics commonly used in PBMC and tumor scRNA-seq.
    """
    if "mt" not in adata.var.columns:
        annotate_gene_sets(adata)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    if logger:
        logger.info("Computed QC metrics: total_counts, n_genes_by_counts, pct_counts_mt/ribo/hb, etc.")
    return adata


def suggest_thresholds(
    adata,
    dataset_type: str = "generic",
    logger: Optional[logging.Logger] = None
) -> Dict[str, Any]:
    """
    Provide starting thresholds (NOT universal truth).
    For tumor, keep mt% permissive by default to avoid dropping malignant/stressed populations too early.
    """
    # robust quantiles
    q = adata.obs[["total_counts", "n_genes_by_counts"]].quantile([0.01, 0.99])
    total_lo, total_hi = float(q.loc[0.01, "total_counts"]), float(q.loc[0.99, "total_counts"])
    genes_lo, genes_hi = float(q.loc[0.01, "n_genes_by_counts"]), float(q.loc[0.99, "n_genes_by_counts"])

    if dataset_type == "pbmc":
        mt_max = 10.0  # typical starting point; tune per dataset
    elif dataset_type == "tumor":
        mt_max = 20.0  # more permissive start; refine later after cell-type/malignant separation
    else:
        mt_max = 15.0

    ribo_hi = float(adata.obs["pct_counts_ribo"].quantile(0.99)) if "pct_counts_ribo" in adata.obs else None
    hb_hi = float(adata.obs["pct_counts_hb"].quantile(0.99)) if "pct_counts_hb" in adata.obs else None

    thresholds = {
        "min_counts": max(200.0, total_lo),
        "max_counts": total_hi,
        "min_genes": max(200.0, genes_lo),
        "max_genes": genes_hi,
        "max_pct_mt": mt_max,
        "max_pct_ribo": ribo_hi,
        "max_pct_hb": hb_hi,
    }
    if logger:
        logger.info(f"Suggested QC thresholds (starter): {thresholds}")
    return thresholds


def filter_cells_qc(
    adata,
    min_counts: Optional[float],
    max_counts: Optional[float],
    min_genes: Optional[float],
    max_genes: Optional[float],
    max_pct_mt: Optional[float],
    max_pct_ribo: Optional[float] = None,
    max_pct_hb: Optional[float] = None,
    logger: Optional[logging.Logger] = None
):
    """
    QC-based cell filtering using common metrics.
    """
    before = adata.n_obs
    mask = np.ones(adata.n_obs, dtype=bool)
    if min_counts is not None:
        mask = mask & (adata.obs["total_counts"] >= min_counts)
    if max_counts is not None:
        mask = mask & (adata.obs["total_counts"] <= max_counts)
    if min_genes is not None:
        mask = mask & (adata.obs["n_genes_by_counts"] >= min_genes)
    if max_genes is not None:
        mask = mask & (adata.obs["n_genes_by_counts"] <= max_genes)
    if max_pct_mt is not None:
        mask = mask & (adata.obs["pct_counts_mt"] <= max_pct_mt)
    if max_pct_ribo is not None:
        mask = mask & (adata.obs["pct_counts_ribo"] <= max_pct_ribo)
    if max_pct_hb is not None:
        mask = mask & (adata.obs["pct_counts_hb"] <= max_pct_hb)
    adata = adata[mask].copy()
    after = adata.n_obs
    if logger:
        logger.info(
            "Filtered cells by QC: "
            f"{before} -> {after} (removed {before-after}) | "
            f"min_counts={min_counts}, max_counts={max_counts}, "
            f"min_genes={min_genes}, max_genes={max_genes}, max_pct_mt={max_pct_mt}, "
            f"max_pct_ribo={max_pct_ribo}, max_pct_hb={max_pct_hb}"
        )
    return adata


def basic_qc_summary(adata) -> pd.DataFrame:
    cols = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
    ]
    cols = [c for c in cols if c in adata.obs.columns]
    return adata.obs[cols].describe().T


def per_sample_qc_summary(adata, sample_col: str = "sample") -> pd.DataFrame:
    cols = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
    ]
    cols = [c for c in cols if c in adata.obs.columns]
    if sample_col not in adata.obs.columns:
        raise ValueError(f"sample_col '{sample_col}' not found in adata.obs")
    g = adata.obs.groupby(sample_col)[cols]
    out = g.agg(["median", "mean", "count"])
    return out
