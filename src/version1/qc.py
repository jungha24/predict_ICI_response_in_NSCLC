from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Literal
from pathlib import Path
import logging

import numpy as np
import scanpy as sc

from .utils_io import (
    RunConfig, ensure_dirs, setup_logger, Timer,
    save_checkpoint, assert_basic_sanity
)
from .utils_qc import (
    annotate_gene_sets, compute_qc_metrics,
    suggest_thresholds, filter_cells_qc
)

# Optional dependencies
try:
    import scrublet as scr
    _HAS_SCRUBLET = True
except Exception:
    _HAS_SCRUBLET = False


@dataclass
class QCParams:
    dataset_type: Optional[Literal["pbmc", "tumor", "generic"]] = None
    min_counts: Optional[float] = None
    use_min_counts_filter: Optional[bool] = None
    max_counts: Optional[float] = None
    use_max_counts_filter: Optional[bool] = None
    min_genes: Optional[float] = None
    use_min_genes_filter: Optional[bool] = None
    max_genes: Optional[float] = None
    max_pct_mt: Optional[float] = None
    use_max_pct_mt_filter: Optional[bool] = None
    use_max_genes_filter: Optional[bool] = None

    max_pct_ribo: Optional[float] = None
    use_max_pct_ribo_filter: Optional[bool] = None

    max_pct_hb: Optional[float] = None
    use_max_pct_hb_filter: Optional[Literal["auto", "on", "off"]] = None

    run_doublet: Optional[bool] = None
    expected_doublet_rate: Optional[float] = None

    keep_raw_layer: bool = True


def run_qc(
    adata,
    cfg: RunConfig,
    params: Optional[QCParams] = None,
    logger: Optional[logging.Logger] = None,
    save_name: str = "adata_qc",
):
    """
    QC pipeline:
    - sanity checks
    - qc metrics (mt/ribo/hb)
    - threshold suggestion + filtering
    - optional doublet scoring (Scrublet)
    - checkpoint save
    """
    ensure_dirs(cfg)
    logger = logger or setup_logger(name="qc")
    params = params or QCParams()

    np.random.seed(cfg.seed)
    sc.settings.verbosity = 2
    sc.settings.figdir = str(Path(cfg.project_dir) / cfg.fig_dir)

    assert_basic_sanity(adata, cfg, logger)

    # Preserve raw counts if needed
    if params.keep_raw_layer:
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()

    with Timer(logger, "QC metrics"):
        annotate_gene_sets(adata)
        compute_qc_metrics(adata, logger=logger)

    dataset_type_for_suggest = params.dataset_type if params.dataset_type is not None else cfg.dataset_type
    use_max_genes_filter = (
        params.use_max_genes_filter if params.use_max_genes_filter is not None else cfg.use_max_genes_filter
    )
    use_min_counts_filter = (
        params.use_min_counts_filter if params.use_min_counts_filter is not None else cfg.use_min_counts_filter
    )
    use_max_counts_filter = (
        params.use_max_counts_filter if params.use_max_counts_filter is not None else cfg.use_max_counts_filter
    )
    use_min_genes_filter = (
        params.use_min_genes_filter if params.use_min_genes_filter is not None else cfg.use_min_genes_filter
    )
    use_max_pct_mt_filter = (
        params.use_max_pct_mt_filter if params.use_max_pct_mt_filter is not None else cfg.use_max_pct_mt_filter
    )
    use_max_pct_ribo_filter = (
        params.use_max_pct_ribo_filter if params.use_max_pct_ribo_filter is not None else cfg.use_max_pct_ribo_filter
    )
    hb_mode = params.use_max_pct_hb_filter if params.use_max_pct_hb_filter is not None else cfg.use_max_pct_hb_filter
    if hb_mode not in {"auto", "on", "off"}:
        raise ValueError("use_max_pct_hb_filter must be one of: auto, on, off")
    use_max_pct_hb_filter = hb_mode == "on" or (hb_mode == "auto" and dataset_type_for_suggest == "pbmc")

    min_counts = (params.min_counts if params.min_counts is not None else cfg.min_counts) if use_min_counts_filter else None
    max_counts = (params.max_counts if params.max_counts is not None else cfg.max_counts) if use_max_counts_filter else None
    min_genes = (params.min_genes if params.min_genes is not None else cfg.min_genes) if use_min_genes_filter else None
    max_genes = (params.max_genes if params.max_genes is not None else cfg.max_genes) if use_max_genes_filter else None
    max_pct_mt = (params.max_pct_mt if params.max_pct_mt is not None else cfg.max_pct_mt) if use_max_pct_mt_filter else None
    max_pct_ribo = (
        (params.max_pct_ribo if params.max_pct_ribo is not None else cfg.max_pct_ribo)
        if use_max_pct_ribo_filter else None
    )
    max_pct_hb = (
        (params.max_pct_hb if params.max_pct_hb is not None else cfg.max_pct_hb)
        if use_max_pct_hb_filter else None
    )

    needs_suggest = (
        (use_min_counts_filter and min_counts is None) or
        (use_max_counts_filter and max_counts is None) or
        (use_min_genes_filter and min_genes is None) or
        (use_max_pct_mt_filter and max_pct_mt is None) or
        (use_max_genes_filter and max_genes is None) or
        (use_max_pct_ribo_filter and max_pct_ribo is None) or
        (use_max_pct_hb_filter and max_pct_hb is None)
    )

    if needs_suggest:
        th = suggest_thresholds(adata, dataset_type=dataset_type_for_suggest, logger=logger)
        if use_min_counts_filter:
            min_counts = min_counts if min_counts is not None else th["min_counts"]
        if use_max_counts_filter:
            max_counts = max_counts if max_counts is not None else th["max_counts"]
        if use_min_genes_filter:
            min_genes = min_genes if min_genes is not None else th["min_genes"]
        if use_max_genes_filter:
            max_genes = max_genes if max_genes is not None else th["max_genes"]
        if use_max_pct_mt_filter:
            max_pct_mt = max_pct_mt if max_pct_mt is not None else th["max_pct_mt"]
        if use_max_pct_ribo_filter:
            max_pct_ribo = max_pct_ribo if max_pct_ribo is not None else th.get("max_pct_ribo")
        if use_max_pct_hb_filter:
            max_pct_hb = max_pct_hb if max_pct_hb is not None else th.get("max_pct_hb")
    else:
        logger.info("Using configured QC thresholds (suggest_thresholds skipped).")

    if use_max_pct_ribo_filter and max_pct_ribo is None:
        raise ValueError("max_pct_ribo filter is enabled but no threshold is available. Set max_pct_ribo explicitly.")
    if use_max_pct_hb_filter and max_pct_hb is None:
        raise ValueError("max_pct_hb filter is enabled but no threshold is available. Set max_pct_hb explicitly.")

    with Timer(logger, "QC filter cells"):
        adata = filter_cells_qc(
            adata,
            min_counts=min_counts,
            max_counts=max_counts,
            min_genes=min_genes,
            max_genes=max_genes,
            max_pct_mt=max_pct_mt,
            max_pct_ribo=max_pct_ribo,
            max_pct_hb=max_pct_hb,
            logger=logger,
        )

    run_doublet = params.run_doublet if params.run_doublet is not None else cfg.run_doublet
    if run_doublet:
        if not _HAS_SCRUBLET:
            logger.warning("Scrublet not installed. Skipping doublet detection.")
        else:
            with Timer(logger, "Doublet detection (Scrublet)"):
                X = adata.layers["counts"] if "counts" in adata.layers else adata.X
                expected_rate = (
                    params.expected_doublet_rate
                    if params.expected_doublet_rate is not None
                    else cfg.expected_doublet_rate
                )
                scrub = scr.Scrublet(X, expected_doublet_rate=expected_rate)
                doublet_scores, predicted_doublets = scrub.scrub_doublets()
                adata.obs["doublet_score"] = doublet_scores
                adata.obs["predicted_doublet"] = predicted_doublets
                logger.info(f"Doublets predicted: {predicted_doublets.sum()} / {adata.n_obs}")

    save_checkpoint(adata, cfg, save_name, logger=logger)
    return adata
