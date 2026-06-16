from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


GROUP_SEP = "|||"


def make_group_key(analysis_id: str, celltype: str) -> str:
    return f"{analysis_id}{GROUP_SEP}{celltype}"


def split_group_key(group_key: str) -> Tuple[str, str]:
    analysis_id, celltype = group_key.rsplit(GROUP_SEP, 1)
    return analysis_id, celltype


def first_non_null(series: pd.Series):
    vals = series.dropna()
    if len(vals) == 0:
        return np.nan
    return vals.iloc[0]


def _empty_valid_groups() -> pd.DataFrame:
    return pd.DataFrame({
        "analysis_id": pd.Series(dtype="object"),
        "celltype": pd.Series(dtype="object"),
        "n_cells": pd.Series(dtype="int64"),
        "group_key": pd.Series(dtype="object"),
    })


def _empty_pseudobulk_meta(meta_use: pd.DataFrame) -> pd.DataFrame:
    meta_cols = ["group_key", "analysis_id", "patient_id", "source", "batch", "timepoint", "celltype", "n_cells"]
    if "age" in meta_use.columns:
        meta_cols.append("age")
    if "sex" in meta_use.columns:
        meta_cols.append("sex")
    return pd.DataFrame(columns=meta_cols)


def get_valid_pseudobulk_groups(meta_use: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Only pseudobulk celltypes (PCA/NMF union) are considered for pseudobulk eligibility.
    """
    pb_celltypes = list(dict.fromkeys(
        cfg["celltypes"].get("pseudobulk_pca", []) +
        cfg["celltypes"].get("pseudobulk_nmf", [])
    ))
    if len(pb_celltypes) == 0:
        return _empty_valid_groups()

    min_cells = int(cfg["pseudobulk"]["min_cells_per_patient_celltype"])

    meta_pb = meta_use.loc[meta_use["celltype"].isin(pb_celltypes)].copy()
    if meta_pb.empty:
        return _empty_valid_groups()

    group_counts = (
        meta_pb.groupby(["analysis_id", "celltype"])
        .size()
        .reset_index(name="n_cells")
    )
    group_counts = group_counts.loc[group_counts["n_cells"] >= min_cells].copy()
    if group_counts.empty:
        return _empty_valid_groups()

    group_counts["group_key"] = group_counts.apply(
        lambda r: make_group_key(r["analysis_id"], r["celltype"]), axis=1
    )
    return group_counts


def aggregate_pseudobulk(
    counts_use: sp.csr_matrix,
    meta_use: pd.DataFrame,
    valid_groups: pd.DataFrame,
) -> Tuple[sp.csr_matrix, List[str], pd.DataFrame]:
    """
    Aggregate cells into genes x pseudobulk-groups sparse matrix.

    counts_use: genes x cells
    meta_use: cell-level metadata aligned to counts_use columns
    valid_groups: dataframe from get_valid_pseudobulk_groups()

    returns:
      pb_counts: genes x pseudobulk-groups
      pb_group_names: list of group keys
      pb_meta: group-level metadata aligned to pseudobulk columns
    """
    if valid_groups.empty:
        return sp.csr_matrix((counts_use.shape[0], 0), dtype=counts_use.dtype), [], _empty_pseudobulk_meta(meta_use)

    valid_groups = valid_groups.copy()
    valid_groups["group_key"] = valid_groups["group_key"].astype(str)
    valid_group_keys = set(valid_groups["group_key"].tolist())

    meta_pb = meta_use.copy().reset_index(drop=True)
    meta_pb["group_key"] = meta_pb.apply(lambda r: make_group_key(r["analysis_id"], r["celltype"]), axis=1)
    meta_pb = meta_pb.loc[meta_pb["group_key"].isin(valid_group_keys)].copy()
    if meta_pb.empty:
        return sp.csr_matrix((counts_use.shape[0], 0), dtype=counts_use.dtype), [], _empty_pseudobulk_meta(meta_use)

    keep_idx = meta_pb.index.to_numpy()
    counts_pb_input = counts_use[:, keep_idx]

    groups = meta_pb["group_key"].astype("category")
    group_codes = groups.cat.codes.to_numpy()
    n_groups = groups.cat.categories.size

    indicator = sp.csr_matrix(
        (np.ones(len(group_codes)), (np.arange(len(group_codes)), group_codes)),
        shape=(len(group_codes), n_groups),
    )

    pb_counts = counts_pb_input @ indicator
    pb_group_names = groups.cat.categories.tolist()

    meta_cols = ["group_key", "analysis_id", "patient_id", "source", "batch", "timepoint", "celltype"]
    if "age" in meta_pb.columns:
        meta_cols.append("age")
    if "sex" in meta_pb.columns:
        meta_cols.append("sex")

    pb_meta = (
        meta_pb[meta_cols]
        .groupby("group_key", as_index=False)
        .agg({
            "analysis_id": first_non_null,
            "patient_id": first_non_null,
            "source": first_non_null,
            "batch": first_non_null,
            "timepoint": first_non_null,
            "celltype": first_non_null,
            **({"age": first_non_null} if "age" in meta_cols else {}),
            **({"sex": first_non_null} if "sex" in meta_cols else {}),
        })
    )
    pb_meta = pb_meta.merge(valid_groups[["group_key", "n_cells"]], on="group_key", how="left")
    pb_meta = pb_meta.set_index("group_key").reindex(pb_group_names).reset_index()

    return pb_counts.tocsr(), pb_group_names, pb_meta
