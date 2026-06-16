from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


CANONICAL_META_KEYS = ["patient_id", "source", "batch", "timepoint", "celltype", "age", "sex"]


def _get_metadata_mapping(cfg: dict) -> Dict[str, Optional[str]]:
    meta_cfg = cfg.get("metadata_cols", {})
    return {
        "patient_id": meta_cfg.get("patient_id"),
        "source": meta_cfg.get("source"),
        "batch": meta_cfg.get("batch"),
        "timepoint": meta_cfg.get("timepoint"),
        "celltype": meta_cfg.get("celltype"),
        "age": meta_cfg.get("age"),
        "sex": meta_cfg.get("sex"),
    }


def canonicalize_metadata(meta: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Rename user-provided metadata column names to internal canonical names:
    patient_id, source, batch, timepoint, celltype, age, sex
    """
    meta = meta.copy()
    mapping = _get_metadata_mapping(cfg)

    if "cell_id" not in meta.columns:
        raise ValueError("cell_meta.parquet must contain a 'cell_id' column.")

    rename_dict = {}
    for canon_key, original_col in mapping.items():
        if original_col is None:
            continue
        if original_col not in meta.columns:
            raise ValueError(
                f"Configured metadata column '{original_col}' for '{canon_key}' was not found in cell_meta.parquet."
            )
        rename_dict[original_col] = canon_key

    meta = meta.rename(columns=rename_dict)

    required = ["cell_id", "patient_id", "source", "batch", "timepoint", "celltype"]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"Missing required canonical metadata columns after rename: {missing}")

    meta["cell_id"] = meta["cell_id"].astype(str)
    meta["patient_id"] = meta["patient_id"].astype(str)
    meta["source"] = meta["source"].astype(str)
    meta["batch"] = meta["batch"].astype(str)
    meta["timepoint"] = meta["timepoint"].astype(str)
    meta["celltype"] = meta["celltype"].astype(str)

    if "age" in meta.columns:
        meta["age"] = pd.to_numeric(meta["age"], errors="coerce")

    if "sex" in meta.columns:
        meta["sex"] = meta["sex"].astype(str)

    return meta


def read_selected_patients(file_path: str | Path) -> List[str]:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Selected patient list not found: {file_path}")

    with open(file_path, "r", encoding="utf-8") as f:
        patient_ids = [line.strip() for line in f if line.strip()]

    patient_ids = list(dict.fromkeys(patient_ids))
    return patient_ids


def get_required_celltypes(cfg: dict) -> List[str]:
    ct_cfg = cfg.get("celltypes", {})
    union_ct = []
    for key in ["proportion", "pseudobulk_pca", "pseudobulk_nmf"]:
        vals = ct_cfg.get(key, [])
        union_ct.extend(vals)
    return list(dict.fromkeys(union_ct))


def _sample_healthy_ids(meta: pd.DataFrame, cfg: dict) -> List[str]:
    hs_cfg = cfg.get("healthy_sampling", {})
    if not hs_cfg.get("enabled", True):
        return []

    healthy_value = cfg["source_values"]["healthy"]
    n_healthy = int(hs_cfg.get("n_healthy", 0))
    random_seed = int(hs_cfg.get("random_seed", cfg.get("seed", 123)))

    healthy_ids = (
        meta.loc[meta["source"] == healthy_value, "patient_id"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    rng = np.random.default_rng(random_seed)
    if len(healthy_ids) == 0 or n_healthy <= 0:
        return []

    if len(healthy_ids) <= n_healthy:
        return sorted(healthy_ids)

    sampled = rng.choice(healthy_ids, size=n_healthy, replace=False)
    return sorted(sampled.tolist())


def build_discovery_cohort(
    counts,
    meta: pd.DataFrame,
    cfg: dict,
) -> Tuple[object, pd.DataFrame, dict]:
    """
    Build subtype discovery cohort using:
      - disease == configured disease source
      - timepoint == baseline value
      - optionally restricted to selected patient list
      - plus randomly sampled healthy donors
      - only celltypes required for downstream steps
    """
    meta = canonicalize_metadata(meta, cfg)

    disease_value = cfg["source_values"]["disease"]
    healthy_value = cfg["source_values"]["healthy"]
    baseline_value = cfg["baseline"]["value"]

    required_celltypes = set(get_required_celltypes(cfg))

    # patient selection
    ps_cfg = cfg.get("patient_selection", {})
    selected_patients: Optional[List[str]] = None
    if ps_cfg.get("enabled", False):
        selected_patients = read_selected_patients(ps_cfg["file"])

    healthy_ids = _sample_healthy_ids(meta, cfg)

    disease_mask = (
        (meta["source"] == disease_value)
        & (meta["timepoint"] == baseline_value)
    )

    if selected_patients is not None and ps_cfg.get("apply_to_disease_only", True):
        disease_mask &= meta["patient_id"].isin(selected_patients)

    healthy_mask = pd.Series(False, index=meta.index)
    if len(healthy_ids) > 0:
        healthy_mask = (
            (meta["source"] == healthy_value)
            & (meta["patient_id"].isin(healthy_ids))
        )

    celltype_mask = meta["celltype"].isin(required_celltypes)

    keep_mask = (disease_mask | healthy_mask) & celltype_mask

    meta_use = meta.loc[keep_mask].copy().reset_index(drop=True)
    keep_idx = np.where(keep_mask.values)[0]
    counts_use = counts[:, keep_idx]

    meta_use["analysis_id"] = meta_use["source"].astype(str) + "__" + meta_use["patient_id"].astype(str)

    cohort_info = {
        "selected_patients_file": ps_cfg.get("file"),
        "selected_patients_n": 0 if selected_patients is None else len(selected_patients),
        "healthy_sampled_ids": healthy_ids,
        "n_cells_total": int(meta_use.shape[0]),
        "n_patients_total": int(meta_use["analysis_id"].nunique()),
        "n_disease_patients": int(meta_use.loc[meta_use["source"] == disease_value, "analysis_id"].nunique()),
        "n_healthy_patients": int(meta_use.loc[meta_use["source"] == healthy_value, "analysis_id"].nunique()),
        "required_celltypes": sorted(required_celltypes),
    }

    return counts_use, meta_use, cohort_info


def summarize_discovery_cohort(meta_use: pd.DataFrame) -> dict:
    out = {
        "patients_by_source": (
            meta_use.groupby("source")["analysis_id"]
            .nunique()
            .to_dict()
        ),
        "cells_by_source": (
            meta_use.groupby("source")
            .size()
            .to_dict()
        ),
        "cells_by_celltype": (
            meta_use.groupby("celltype")
            .size()
            .sort_values(ascending=False)
            .to_dict()
        ),
        "patient_celltype_table": (
            meta_use.groupby(["analysis_id", "celltype"])
            .size()
            .reset_index(name="n_cells")
            .head(20)
            .to_dict(orient="records")
        ),
    }
    return out
