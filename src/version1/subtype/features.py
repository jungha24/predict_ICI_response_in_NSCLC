from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.linalg import helmert
from sklearn.decomposition import PCA, NMF
from sklearn.preprocessing import StandardScaler

from .pseudobulk import split_group_key


def log_cpm(mat: sp.spmatrix) -> np.ndarray:
    """
    mat: genes x samples sparse matrix
    returns: genes x samples dense log1p(CPM)
    """
    libsize = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    libsize[libsize == 0] = 1.0
    scale = 1e6 / libsize
    cpm = mat @ sp.diags(scale)
    return np.log1p(cpm.toarray())


def nonlog_cpm(mat: sp.spmatrix) -> np.ndarray:
    """
    mat: genes x samples sparse matrix
    returns: genes x samples dense CPM
    """
    libsize = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    libsize[libsize == 0] = 1.0
    scale = 1e6 / libsize
    cpm = mat @ sp.diags(scale)
    return cpm.toarray()


def first_non_null(series: pd.Series):
    vals = series.dropna()
    if len(vals) == 0:
        return np.nan
    return vals.iloc[0]


def _resolve_proportion_cfg(cfg: dict) -> Tuple[bool, str, float]:
    prop_cfg = cfg.get("proportion_features", {})
    enabled = bool(prop_cfg.get("enabled", True))
    transform = str(prop_cfg.get("transform", "clr")).strip().lower()
    pseudocount = float(prop_cfg.get("pseudocount", 0.5))

    if transform not in {"raw", "clr", "ilr"}:
        raise ValueError("proportion_features.transform must be one of: raw, clr, ilr")
    if transform in {"clr", "ilr"} and pseudocount <= 0:
        raise ValueError("proportion_features.pseudocount must be > 0 when using CLR/ILR")

    return enabled, transform, pseudocount


def _prefix_prop_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.reset_index().copy()
    out.columns = ["analysis_id"] + [f"prop_{c}" for c in df.columns]
    return out


def _make_clr_features(prop_df: pd.DataFrame) -> pd.DataFrame:
    log_prop = np.log(prop_df.to_numpy(dtype=float))
    clr = log_prop - log_prop.mean(axis=1, keepdims=True)
    return pd.DataFrame(clr, index=prop_df.index, columns=prop_df.columns)


def _make_ilr_features(prop_df: pd.DataFrame) -> pd.DataFrame:
    n_parts = prop_df.shape[1]
    if n_parts < 2:
        return pd.DataFrame(index=prop_df.index)

    basis = helmert(n_parts, full=False).T
    clr = _make_clr_features(prop_df).to_numpy(dtype=float)
    ilr = clr @ basis
    columns = [f"ilr{i+1}" for i in range(ilr.shape[1])]
    return pd.DataFrame(ilr, index=prop_df.index, columns=columns)


def build_patient_metadata(meta_use: pd.DataFrame) -> pd.DataFrame:
    keep_cols = ["analysis_id", "patient_id", "source", "batch", "timepoint"]
    if "age" in meta_use.columns:
        keep_cols.append("age")
    if "sex" in meta_use.columns:
        keep_cols.append("sex")

    patient_meta = (
        meta_use[keep_cols]
        .groupby("analysis_id", as_index=False)
        .agg({
            "patient_id": first_non_null,
            "source": first_non_null,
            "batch": first_non_null,
            "timepoint": first_non_null,
            **({"age": first_non_null} if "age" in keep_cols else {}),
            **({"sex": first_non_null} if "sex" in keep_cols else {}),
        })
    )
    return patient_meta


def make_proportion_features(meta_use: pd.DataFrame, cfg: dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Proportions are computed only among configured proportion celltypes.
    """
    prop_celltypes = cfg["celltypes"]["proportion"]
    enabled, transform, pseudocount = _resolve_proportion_cfg(cfg)
    analysis_ids = pd.Index(meta_use["analysis_id"].drop_duplicates(), name="analysis_id")

    if not enabled:
        empty_wide = pd.DataFrame({"analysis_id": analysis_ids.to_numpy()})
        empty_long = pd.DataFrame(columns=["analysis_id", "celltype", "n_cells", "prop"])
        return empty_wide.copy(), empty_long, empty_wide.copy()

    meta_prop = meta_use.loc[meta_use["celltype"].isin(prop_celltypes)].copy()

    prop_long = (
        meta_prop.groupby(["analysis_id", "celltype"])
        .size()
        .reset_index(name="n_cells")
    )
    prop_long["prop"] = prop_long.groupby("analysis_id")["n_cells"].transform(lambda x: x / x.sum())

    count_wide = (
        prop_long.pivot(index="analysis_id", columns="celltype", values="n_cells")
        .reindex(index=analysis_ids, columns=prop_celltypes, fill_value=0)
    )
    raw_prop_wide = count_wide.div(count_wide.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    diagnostic_wide = _prefix_prop_columns(raw_prop_wide)

    if transform == "raw":
        transformed = raw_prop_wide
        diagnostic_wide = _prefix_prop_columns(transformed)
    else:
        smoothed_prop = (count_wide + pseudocount).div((count_wide + pseudocount).sum(axis=1), axis=0)
        if transform == "clr":
            transformed = _make_clr_features(smoothed_prop)
            diagnostic_wide = _prefix_prop_columns(transformed)
        else:
            transformed = _make_ilr_features(smoothed_prop)

    prop_wide = _prefix_prop_columns(transformed)

    return prop_wide, prop_long, diagnostic_wide


def make_pseudobulk_pca_features(
    pb_counts: sp.spmatrix,
    pb_group_names: List[str],
    genes: List[str],
    cfg: dict,
    expr_override: Optional[np.ndarray] = None,
) -> Optional[pd.DataFrame]:
    """
    pb_counts: genes x groups
    """
    pca_celltypes = cfg["celltypes"]["pseudobulk_pca"]
    npcs = int(cfg["pca_features"]["npcs_per_celltype"])
    min_patients = int(cfg["pca_features"]["min_patients_per_celltype_pca"])
    expr_all = expr_override if expr_override is not None else log_cpm(pb_counts)

    out = []

    for ct in pca_celltypes:
        cols = [i for i, g in enumerate(pb_group_names) if split_group_key(g)[1] == ct]
        if len(cols) < min_patients:
            continue

        expr = expr_all[:, cols]

        # simple gene filtering: keep genes with non-zero variance
        gene_var = expr.var(axis=1)
        keep = gene_var > 0
        expr = expr[keep, :]
        if expr.shape[0] < 50:
            continue

        sample_ids = [split_group_key(pb_group_names[i])[0] for i in cols]
        x = expr.T  # patients x genes

        x_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(x)
        n_comp = min(npcs, x_scaled.shape[0], x_scaled.shape[1])
        if n_comp < 1:
            continue

        pca = PCA(n_components=n_comp, random_state=cfg.get("seed", 123))
        pcs = pca.fit_transform(x_scaled)

        df = pd.DataFrame(
            pcs,
            columns=[f"{ct}_PC{i+1}" for i in range(pcs.shape[1])]
        )
        df["analysis_id"] = sample_ids
        out.append(df)

    if len(out) == 0:
        return None

    merged = out[0]
    for df in out[1:]:
        merged = merged.merge(df, on="analysis_id", how="outer")
    return merged


def make_pseudobulk_nmf_features(
    pb_counts: sp.spmatrix,
    pb_group_names: List[str],
    genes: List[str],
    cfg: dict,
    log_expr_override: Optional[np.ndarray] = None,
    nmf_expr_override: Optional[np.ndarray] = None,
) -> Tuple[Optional[pd.DataFrame], Dict[str, dict]]:
    """
    NMF is run only if enabled == true.
    Returns:
      - merged activity dataframe
      - per-celltype top genes metadata
    """
    nmf_cfg = cfg.get("nmf_features", {})
    if not nmf_cfg.get("enabled", True):
        return None, {}

    nmf_celltypes = cfg["celltypes"]["pseudobulk_nmf"]
    rank = int(nmf_cfg["rank"])
    min_patients = int(nmf_cfg["min_patients_per_celltype_nmf"])
    min_genes = int(nmf_cfg["min_genes_nmf"])
    top_var_genes = int(nmf_cfg["nmf_top_var_genes"])

    out = []
    nmf_meta = {}

    gene_names = np.array(genes)
    expr_log_all = log_expr_override if log_expr_override is not None else log_cpm(pb_counts)
    expr_nmf_all = nmf_expr_override if nmf_expr_override is not None else nonlog_cpm(pb_counts)

    for ct in nmf_celltypes:
        cols = [i for i, g in enumerate(pb_group_names) if split_group_key(g)[1] == ct]
        if len(cols) < min_patients:
            continue

        expr_log = expr_log_all[:, cols]
        gene_var = expr_log.var(axis=1)
        keep_var = np.argsort(-gene_var)[: min(top_var_genes, len(gene_var))]
        expr_nmf = np.clip(expr_nmf_all[np.ix_(keep_var, cols)], a_min=0.0, a_max=None) + 1e-8  # genes x patients

        if expr_nmf.shape[0] < min_genes:
            continue

        x = expr_nmf.T  # patients x genes, non-negative
        n_comp = min(rank, x.shape[0], x.shape[1])
        if n_comp < 1:
            continue

        model = NMF(
            n_components=n_comp,
            init="nndsvda",
            random_state=cfg.get("seed", 123),
            max_iter=1000,
        )
        W = model.fit_transform(x)    # patients x modules
        H = model.components_         # modules x genes

        sample_ids = [split_group_key(pb_group_names[i])[0] for i in cols]

        df = pd.DataFrame(
            W,
            columns=[f"{ct}_NMF{i+1}" for i in range(W.shape[1])]
        )
        df["analysis_id"] = sample_ids
        out.append(df)

        kept_gene_names = gene_names[keep_var]
        top_genes = {}
        for i in range(H.shape[0]):
            order = np.argsort(-H[i])[:30]
            top_genes[f"{ct}_NMF{i+1}"] = kept_gene_names[order].tolist()

        nmf_meta[ct] = {
            "rank_used": int(n_comp),
            "top_genes": top_genes,
        }

    if len(out) == 0:
        return None, nmf_meta

    merged = out[0]
    for df in out[1:]:
        merged = merged.merge(df, on="analysis_id", how="outer")

    return merged, nmf_meta


def fill_missing_values(patient_features: pd.DataFrame) -> pd.DataFrame:
    """
    - proportion features -> fill 0
    - transcript features -> fill median
    """
    out = patient_features.copy()
    meta_like = {"analysis_id", "patient_id", "source", "batch", "timepoint", "age", "sex"}

    for col in out.columns:
        if col in meta_like:
            continue

        if col.startswith("prop_"):
            out[col] = out[col].fillna(0.0)
        else:
            med = out[col].median(skipna=True)
            if pd.isna(med):
                med = 0.0
            out[col] = out[col].fillna(med)

    return out
