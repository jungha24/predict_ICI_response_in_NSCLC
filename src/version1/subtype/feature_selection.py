from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


META_EXCLUDE = {"analysis_id", "patient_id", "source", "batch", "timepoint", "age", "sex"}


def get_numeric_feature_columns(df: pd.DataFrame) -> List[str]:
    cols = [c for c in df.columns if c not in META_EXCLUDE and pd.api.types.is_numeric_dtype(df[c])]
    return cols


def extract_celltype_name(feature: str) -> str | None:
    if feature.startswith("prop_"):
        return feature.replace("prop_", "", 1)
    if "_PC" in feature:
        return feature.split("_PC")[0]
    if "_NMF" in feature:
        return feature.split("_NMF")[0]
    return None


def infer_feature_type(feature: str) -> str | None:
    if feature.startswith("prop_"):
        return "proportion"
    if "_PC" in feature:
        return "pca"
    if "_NMF" in feature:
        return "nmf"
    return None


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) != len(y):
        return np.nan
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    x_use = x[mask]
    y_use = y[mask]
    if np.nanstd(x_use) == 0 or np.nanstd(y_use) == 0:
        return np.nan
    return float(np.corrcoef(x_use, y_use)[0, 1])


def _complete_case_xy(df: pd.DataFrame, x_col: str, y_col: str) -> tuple[np.ndarray, np.ndarray]:
    x = pd.to_numeric(df[x_col], errors="coerce")
    y = pd.to_numeric(df[y_col], errors="coerce")
    mask = x.notna() & y.notna()
    return x.loc[mask].to_numpy().reshape(-1, 1), y.loc[mask].to_numpy()


def _normalize_anchor_map(anchor_map: Any) -> Dict[str, dict]:
    if anchor_map is None:
        return {}
    if isinstance(anchor_map, dict):
        return {str(k): dict(v or {}) for k, v in anchor_map.items()}
    raise ValueError("feature_selection.anchor_feature_map must be a mapping from anchor feature to config")


def _resolve_configured_feature_celltypes(cfg: dict | None) -> Dict[str, set[str]]:
    celltype_cfg = (cfg or {}).get("celltypes", {})
    return {
        "pca": {str(x) for x in celltype_cfg.get("pseudobulk_pca", [])},
        "nmf": {str(x) for x in celltype_cfg.get("pseudobulk_nmf", [])},
    }


def _resolve_anchor_specs(df: pd.DataFrame, cfg: dict | None = None) -> List[dict]:
    feat_cols = get_numeric_feature_columns(df)
    fs_cfg = (cfg or {}).get("feature_selection", {})
    anchor_map = _normalize_anchor_map(fs_cfg.get("anchor_feature_map"))
    configured_feature_celltypes = _resolve_configured_feature_celltypes(cfg)

    specs: List[dict] = []
    if anchor_map:
        for anchor_feature, spec in anchor_map.items():
            if anchor_feature not in feat_cols:
                continue
            related_celltypes = [str(x) for x in spec.get("related_celltypes", [])]
            feature_types = {str(x).lower() for x in spec.get("feature_types", ["pca", "nmf"])}
            if not related_celltypes:
                continue
            allowed_celltypes: set[str] = set()
            for feature_type in feature_types:
                allowed_celltypes.update(configured_feature_celltypes.get(feature_type, set()))
            if allowed_celltypes:
                related_celltypes = [celltype for celltype in related_celltypes if celltype in allowed_celltypes]
            if not related_celltypes:
                continue
            specs.append(
                {
                    "anchor_feature": anchor_feature,
                    "anchor_type": str(spec.get("anchor_type", "anchor")),
                    "related_celltypes": related_celltypes,
                    "feature_types": feature_types,
                }
            )
        return specs

    prop_cols = [c for c in feat_cols if c.startswith("prop_")]
    for pcol in prop_cols:
        ct = pcol.replace("prop_", "", 1)
        specs.append(
            {
                "anchor_feature": pcol,
                "anchor_type": "proportion",
                "related_celltypes": [ct],
                "feature_types": {"pca", "nmf"},
            }
        )
    return specs


def _match_related_features(feature_columns: List[str], spec: dict) -> List[str]:
    matched = []
    for feature in feature_columns:
        if feature == spec["anchor_feature"]:
            continue
        feature_type = infer_feature_type(feature)
        celltype = extract_celltype_name(feature)
        if feature_type not in spec["feature_types"]:
            continue
        if celltype not in spec["related_celltypes"]:
            continue
        matched.append(feature)
    return matched


def check_prop_related_features(df: pd.DataFrame, cfg: dict | None = None) -> pd.DataFrame:
    feat_cols = get_numeric_feature_columns(df)
    rows = []
    for spec in _resolve_anchor_specs(df, cfg):
        anchor_feature = spec["anchor_feature"]
        matched = _match_related_features(feat_cols, spec)
        for feature in matched:
            r = safe_corr(pd.to_numeric(df[anchor_feature], errors="coerce").to_numpy(), pd.to_numeric(df[feature], errors="coerce").to_numpy())
            rows.append(
                {
                    "prop_feature": anchor_feature,
                    "anchor_feature": anchor_feature,
                    "anchor_type": spec["anchor_type"],
                    "feature": feature,
                    "celltype": extract_celltype_name(feature),
                    "related_celltypes": "|".join(spec["related_celltypes"]),
                    "feature_type": infer_feature_type(feature),
                    "correlation": r,
                    "abs_correlation": np.abs(r) if pd.notna(r) else np.nan,
                }
            )

    columns = [
        "prop_feature",
        "anchor_feature",
        "anchor_type",
        "feature",
        "celltype",
        "related_celltypes",
        "feature_type",
        "correlation",
        "abs_correlation",
    ]
    if len(rows) == 0:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows).sort_values("abs_correlation", ascending=False)


def check_feature_residual_signal(df: pd.DataFrame, cfg: dict | None = None) -> pd.DataFrame:
    feat_cols = get_numeric_feature_columns(df)
    rows = []
    for spec in _resolve_anchor_specs(df, cfg):
        anchor_feature = spec["anchor_feature"]
        matched = _match_related_features(feat_cols, spec)
        for feature in matched:
            x, y = _complete_case_xy(df, anchor_feature, feature)
            if len(y) < 3 or np.nanstd(y) == 0 or np.nanstd(x[:, 0]) == 0:
                rows.append(
                    {
                        "prop_feature": anchor_feature,
                        "anchor_feature": anchor_feature,
                        "anchor_type": spec["anchor_type"],
                        "feature": feature,
                        "celltype": extract_celltype_name(feature),
                        "related_celltypes": "|".join(spec["related_celltypes"]),
                        "feature_type": infer_feature_type(feature),
                        "r_squared": np.nan,
                        "raw_sd": float(np.nanstd(y)) if len(y) else np.nan,
                        "residual_sd": np.nan,
                        "residual_ratio": np.nan,
                    }
                )
                continue

            model = LinearRegression().fit(x, y)
            pred = model.predict(x)
            resid = y - pred
            raw_sd = float(np.std(y))
            residual_sd = float(np.std(resid))
            rows.append(
                {
                    "prop_feature": anchor_feature,
                    "anchor_feature": anchor_feature,
                    "anchor_type": spec["anchor_type"],
                    "feature": feature,
                    "celltype": extract_celltype_name(feature),
                    "related_celltypes": "|".join(spec["related_celltypes"]),
                    "feature_type": infer_feature_type(feature),
                    "r_squared": float(model.score(x, y)),
                    "raw_sd": raw_sd,
                    "residual_sd": residual_sd,
                    "residual_ratio": (residual_sd / raw_sd) if raw_sd > 0 else np.nan,
                }
            )

    columns = [
        "prop_feature",
        "anchor_feature",
        "anchor_type",
        "feature",
        "celltype",
        "related_celltypes",
        "feature_type",
        "r_squared",
        "raw_sd",
        "residual_sd",
        "residual_ratio",
    ]
    if len(rows) == 0:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows).sort_values("r_squared", ascending=False)


def check_pca_nmf_redundancy(df: pd.DataFrame) -> pd.DataFrame:
    feat_cols = get_numeric_feature_columns(df)
    pca_cols = [c for c in feat_cols if "_PC" in c]
    nmf_cols = [c for c in feat_cols if "_NMF" in c]

    rows = []
    for pcol in pca_cols:
        ct = extract_celltype_name(pcol)
        matched = [n for n in nmf_cols if extract_celltype_name(n) == ct]
        for ncol in matched:
            r = safe_corr(pd.to_numeric(df[pcol], errors="coerce").to_numpy(), pd.to_numeric(df[ncol], errors="coerce").to_numpy())
            rows.append({
                "celltype": ct,
                "pca_feature": pcol,
                "nmf_feature": ncol,
                "correlation": r,
                "abs_correlation": np.abs(r) if pd.notna(r) else np.nan,
            })

    if len(rows) == 0:
        return pd.DataFrame(columns=["celltype", "pca_feature", "nmf_feature", "correlation", "abs_correlation"])

    return pd.DataFrame(rows).sort_values("abs_correlation", ascending=False)


def select_features_for_clustering(
    patient_features: pd.DataFrame,
    resid_signal_df: pd.DataFrame,
    pca_nmf_red_df: pd.DataFrame,
    cfg: dict,
) -> Tuple[List[str], List[str]]:
    fs_cfg = cfg["feature_selection"]
    high_corr_cutoff = float(fs_cfg["high_corr_cutoff"])
    prop_r2_cutoff = float(fs_cfg["prop_r2_cutoff"])
    residual_ratio_cutoff = float(fs_cfg["residual_ratio_cutoff"])
    prefer = fs_cfg.get("prefer", "pca")
    always_keep_proportions = bool(fs_cfg.get("always_keep_proportions", True))

    all_features = get_numeric_feature_columns(patient_features)
    prop_feats = [c for c in all_features if c.startswith("prop_")]

    drop_feats = set()

    if not resid_signal_df.empty:
        flagged = resid_signal_df[
            (resid_signal_df["r_squared"] >= prop_r2_cutoff) |
            (resid_signal_df["residual_ratio"] <= residual_ratio_cutoff)
        ]
        drop_feats.update(flagged["feature"].tolist())

    if not pca_nmf_red_df.empty:
        redundant = pca_nmf_red_df[pca_nmf_red_df["abs_correlation"] >= high_corr_cutoff]
        for _, row in redundant.iterrows():
            if prefer == "pca":
                drop_feats.add(row["nmf_feature"])
            else:
                drop_feats.add(row["pca_feature"])

    if always_keep_proportions:
        drop_feats -= set(prop_feats)

    selected = [f for f in all_features if f not in drop_feats]
    dropped = sorted([f for f in all_features if f in drop_feats])

    return selected, dropped
