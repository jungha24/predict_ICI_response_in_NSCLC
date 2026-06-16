from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


def resolve_covariates(cfg: dict) -> List[str]:
    """
    Supports both:
      batch_adjustment:
        batch_col: collection_batch, age, sex

    and:
      batch_adjustment:
        covariates: [collection_batch, age, sex]
    """
    ba_cfg = cfg.get("batch_adjustment", {})
    if "covariates" in ba_cfg:
        raw = ba_cfg["covariates"]
    else:
        raw = ba_cfg.get("batch_col", "batch")

    if isinstance(raw, str):
        covariates = [x.strip() for x in raw.split(",") if x.strip()]
    elif isinstance(raw, list):
        covariates = [str(x).strip() for x in raw if str(x).strip()]
    else:
        raise ValueError("batch_adjustment.covariates / batch_col must be str or list")

    meta_cfg = cfg.get("metadata_cols", {})
    original_to_canonical = {}
    for canon_name in ["patient_id", "source", "batch", "timepoint", "celltype", "age", "sex"]:
        original = meta_cfg.get(canon_name)
        if original is not None:
            original_to_canonical[original] = canon_name
        original_to_canonical[canon_name] = canon_name

    resolved = [original_to_canonical.get(c, c) for c in covariates]
    resolved = list(dict.fromkeys(resolved))
    return resolved


def make_design_matrix(df: pd.DataFrame, covariates: List[str]) -> pd.DataFrame:
    design_parts = []

    for cov in covariates:
        if cov not in df.columns:
            continue

        series = df[cov]

        if pd.api.types.is_numeric_dtype(series):
            x = pd.to_numeric(series, errors="coerce")
            x = x.fillna(x.median())
            design_parts.append(pd.DataFrame({cov: x}))
        else:
            x = series.fillna("NA").astype(str)
            dummies = pd.get_dummies(x, prefix=cov, drop_first=True)
            if dummies.shape[1] == 0:
                continue
            design_parts.append(dummies)

    if len(design_parts) == 0:
        return pd.DataFrame(index=df.index)

    return pd.concat(design_parts, axis=1)


def residualize_features(df: pd.DataFrame, feature_cols: list[str], covariates: list[str]) -> pd.DataFrame:
    out = df.copy()
    design = make_design_matrix(out, covariates)

    if design.shape[1] == 0:
        return out

    x = design.values

    for feat in feature_cols:
        y = out[feat].values
        if np.nanstd(y) == 0:
            continue
        model = LinearRegression().fit(x, y)
        pred = model.predict(x)
        resid = y - pred
        out[feat] = resid + np.nanmean(y)

    return out


def residualize_matrix_rows(
    matrix: np.ndarray,
    sample_meta: pd.DataFrame,
    covariates: list[str],
) -> np.ndarray:
    """
    Residualize a dense genes x samples matrix row-wise against the given covariates.

    The adjusted matrix keeps the per-row mean after removing the linear effect of the
    design matrix built from sample-level metadata.
    """
    design = make_design_matrix(sample_meta, covariates)
    out = np.asarray(matrix, dtype=float).copy()

    if out.ndim != 2:
        raise ValueError("matrix must be 2-dimensional (features x samples)")
    if out.shape[1] != sample_meta.shape[0]:
        raise ValueError("matrix columns must match the number of rows in sample_meta")
    if design.shape[1] == 0:
        return out

    x = design.values

    for i in range(out.shape[0]):
        y = out[i, :]
        if np.nanstd(y) == 0:
            continue
        model = LinearRegression().fit(x, y)
        pred = model.predict(x)
        resid = y - pred
        out[i, :] = resid + np.nanmean(y)

    return out
