from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .design import SearchDataset, resolve_clinical_feature_sets
from .io_utils import normalize_id, read_table

ESSENTIAL_CLINICAL_COLUMNS = {
    "patient_id", "sex", "age", "histology", "smoking", "ecog", "drug",
    "pd_event", "pfs_days", "egfr_status", "io_line", "previous_chemo",
    "previous_target", "pd_l1_tps",
}
OPTIONAL_CLINICAL_COLUMNS = {"death_event", "pfs_event", "recist", "benefit"}


def map_sex(value: object) -> object:
    if pd.isna(value):
        return np.nan
    mapping = {"1": "M", "2": "F", "M": "M", "F": "F", "male": "M", "female": "F"}
    return mapping.get(str(value).strip(), value)


def map_smoking(value: object) -> object:
    if pd.isna(value):
        return np.nan
    mapping = {"0": "Never", "1": "Ex", "2": "Current", "Never": "Never", "Ex": "Ex", "Current": "Current"}
    return mapping.get(str(value).strip(), value)


def map_ox(value: object) -> object:
    if pd.isna(value):
        return np.nan
    mapping = {"O": "Yes", "X": "No", "YES": "Yes", "NO": "No", "1": "Yes", "0": "No"}
    return mapping.get(str(value).strip().upper(), value)


def map_drug_class(value: object) -> object:
    if pd.isna(value):
        return np.nan
    drug_name = str(value).strip().lower()
    pd1_drugs = {
        "pembrolizumab", "nivolumab", "cemiplimab", "dostarlimab", "tislelizumab",
        "toripalimab", "sintilimab", "camrelizumab",
    }
    pdl1_drugs = {"atezolizumab", "durvalumab", "avelumab", "sugemalimab", "envafolimab"}
    if drug_name in pd1_drugs:
        return "PD-1"
    if drug_name in pdl1_drugs:
        return "PD-L1"
    return value


def normalize_recist(value: object) -> object:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().upper()
    mapping = {
        "COMPLETE RESPONSE": "CR",
        "PARTIAL RESPONSE": "PR",
        "STABLE DISEASE": "SD",
        "PROGRESSIVE DISEASE": "PD",
        "NON-CR/NON-PD": "NON_CR_NON_PD",
    }
    return mapping.get(text, text)


def coerce_binary_series(series: pd.Series, positive_values: Optional[list[object]] = None, negative_values: Optional[list[object]] = None) -> pd.Series:
    out = pd.Series(np.nan, index=series.index, dtype=float)
    numeric = pd.to_numeric(series, errors="coerce")
    mask_numeric = numeric.notna()
    if mask_numeric.any():
        out.loc[mask_numeric & (numeric == 1)] = 1.0
        out.loc[mask_numeric & (numeric == 0)] = 0.0

    positive_set = {str(value).strip().upper() for value in (positive_values or [1, "1", "YES", "Y", "TRUE", "BENEFIT", "RESPONDER", "DCB"])}
    negative_set = {str(value).strip().upper() for value in (negative_values or [0, "0", "NO", "N", "FALSE", "NO BENEFIT", "NONBENEFIT", "NON-BENEFIT", "NCB"])}
    text = series.astype("string").str.strip().str.upper()
    out.loc[text.isin(positive_set)] = 1.0
    out.loc[text.isin(negative_set)] = 0.0
    return out


def derive_pfs_event(df: pd.DataFrame, derivation_cfg: dict) -> pd.Series:
    pfs_cfg = derivation_cfg.get("pfs_event", {})
    source_columns = [str(column) for column in pfs_cfg.get("source_columns", ["pd_event", "death_event"])]
    usable = [column for column in source_columns if column in df.columns]
    if not usable:
        return pd.Series(np.nan, index=df.index, dtype=float)
    source_df = pd.concat([coerce_binary_series(df[column]) for column in usable], axis=1)
    return source_df.max(axis=1, skipna=True)


def derive_benefit(df: pd.DataFrame, derivation_cfg: dict) -> pd.Series:
    benefit_cfg = derivation_cfg.get("benefit", {})
    recist_column = str(benefit_cfg.get("recist_column", "recist"))
    time_column = str(benefit_cfg.get("time_column", "pfs_days"))
    if recist_column not in df.columns or time_column not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)

    positive_response_values = {str(value).strip().upper() for value in benefit_cfg.get("positive_response_values", ["CR", "PR"])}
    durable_stable_values = {str(value).strip().upper() for value in benefit_cfg.get("durable_stable_values", ["SD"])}
    min_days = float(benefit_cfg.get("min_days", 183))

    recist = df[recist_column].map(normalize_recist)
    pfs_days = pd.to_numeric(df[time_column], errors="coerce")
    out = pd.Series(np.nan, index=df.index, dtype=float)
    is_positive = recist.astype("string").str.upper().isin(positive_response_values)
    is_durable_sd = recist.astype("string").str.upper().isin(durable_stable_values) & (pfs_days >= min_days)
    known_nonbenefit = recist.notna() & ~(is_positive | is_durable_sd)
    out.loc[is_positive | is_durable_sd] = 1.0
    out.loc[known_nonbenefit & pfs_days.notna()] = 0.0
    return out


def derive_restricted_survival(df: pd.DataFrame, spec: dict) -> tuple[pd.Series, pd.Series]:
    source_time_col = str(spec.get("source_time_col", "pfs_days"))
    source_event_col = str(spec.get("source_event_col", "pfs_event"))
    if source_time_col not in df.columns or source_event_col not in df.columns:
        nan_series = pd.Series(np.nan, index=df.index, dtype=float)
        return nan_series.copy(), nan_series.copy()

    horizon_days = float(spec.get("horizon_days", 183))
    source_time = pd.to_numeric(df[source_time_col], errors="coerce")
    source_event = coerce_binary_series(df[source_event_col])
    restricted_time = source_time.copy()
    restricted_time.loc[source_time.notna() & (source_time > horizon_days)] = horizon_days

    restricted_event = pd.Series(np.nan, index=df.index, dtype=float)
    valid = source_time.notna() & source_event.notna()
    restricted_event.loc[valid] = np.where(
        (source_event.loc[valid] == 1.0) & (source_time.loc[valid] <= horizon_days),
        1.0,
        0.0,
    )
    return restricted_time, restricted_event


def derive_binary_survival_endpoint(df: pd.DataFrame, spec: dict) -> tuple[pd.Series, pd.Series]:
    source_time_col = str(spec.get("source_time_col", "pfs_days"))
    source_event_col = str(spec.get("source_event_col", "pfs_event"))
    if source_time_col not in df.columns or source_event_col not in df.columns:
        nan_series = pd.Series(np.nan, index=df.index, dtype=float)
        label_series = pd.Series(pd.NA, index=df.index, dtype="string")
        return nan_series, label_series

    horizon_days = float(spec.get("horizon_days", 183))
    source_time = pd.to_numeric(df[source_time_col], errors="coerce")
    source_event = coerce_binary_series(df[source_event_col])

    out = pd.Series(np.nan, index=df.index, dtype=float)
    labels = pd.Series(pd.NA, index=df.index, dtype="string")

    valid = source_time.notna() & source_event.notna()
    event_by_horizon = valid & (source_event == 1.0) & (source_time <= horizon_days)
    event_free_at_horizon = valid & (
        ((source_event == 0.0) & (source_time >= horizon_days))
        | ((source_event == 1.0) & (source_time > horizon_days))
    )
    censored_before_horizon = valid & (source_event == 0.0) & (source_time < horizon_days)

    positive_label = str(spec.get("positive_label", f"event_by_{int(horizon_days)}d"))
    negative_label = str(spec.get("negative_label", f"event_free_at_{int(horizon_days)}d"))
    censored_label = str(spec.get("censored_label", f"censored_before_{int(horizon_days)}d"))

    out.loc[event_by_horizon] = 1.0
    out.loc[event_free_at_horizon] = 0.0
    labels.loc[event_by_horizon] = positive_label
    labels.loc[event_free_at_horizon] = negative_label
    labels.loc[censored_before_horizon] = censored_label
    return out, labels


def prepare_clinical_df(raw: pd.DataFrame, column_map: dict[str, Optional[str]], derivation_cfg: Optional[dict] = None) -> pd.DataFrame:
    derivation_cfg = derivation_cfg or {}
    missing_required = []
    for std_name in ESSENTIAL_CLINICAL_COLUMNS:
        raw_name = column_map.get(std_name)
        if raw_name is None or raw_name not in raw.columns:
            missing_required.append(std_name if raw_name is None else raw_name)
    if missing_required:
        raise ValueError(f"Missing required clinical columns: {missing_required}")

    df_dict: dict[str, pd.Series] = {}
    for std_name, raw_name in column_map.items():
        if raw_name is None:
            df_dict[std_name] = pd.Series(np.nan, index=raw.index)
        elif raw_name in raw.columns:
            df_dict[std_name] = raw[raw_name]
        elif std_name in OPTIONAL_CLINICAL_COLUMNS:
            df_dict[std_name] = pd.Series(np.nan, index=raw.index)
        else:
            raise ValueError(f"Missing required clinical column: {raw_name}")

    for std_name in OPTIONAL_CLINICAL_COLUMNS:
        if std_name not in df_dict:
            df_dict[std_name] = pd.Series(np.nan, index=raw.index)

    df = pd.DataFrame(df_dict)
    df["patient_id"] = df["patient_id"].map(normalize_id)
    df["sex"] = df["sex"].map(map_sex)
    df["smoking"] = df["smoking"].map(map_smoking)
    df["previous_chemo"] = df["previous_chemo"].map(map_ox)
    df["previous_target"] = df["previous_target"].map(map_ox)
    df["drug"] = df["drug"].astype(str)
    df["drug_class"] = df["drug"].map(map_drug_class)
    df["recist"] = df["recist"].map(normalize_recist)

    numeric_cols = ["age", "ecog", "pd_event", "death_event", "pfs_days", "pfs_event"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if df["pfs_event"].isna().all() or derivation_cfg.get("pfs_event", {}).get("overwrite", True):
        derived_pfs_event = derive_pfs_event(df, derivation_cfg)
        if df["pfs_event"].isna().all():
            df["pfs_event"] = derived_pfs_event
        else:
            df["pfs_event"] = df["pfs_event"].fillna(derived_pfs_event)

    benefit_cfg = derivation_cfg.get("benefit", {})
    benefit_existing = coerce_binary_series(df["benefit"], positive_values=benefit_cfg.get("positive_values"), negative_values=benefit_cfg.get("negative_values"))
    prefer_existing = bool(benefit_cfg.get("prefer_existing", True))
    derived_benefit = derive_benefit(df, derivation_cfg)
    df["benefit"] = benefit_existing.fillna(derived_benefit) if prefer_existing else derived_benefit.fillna(benefit_existing)

    df["pd_event"] = coerce_binary_series(df["pd_event"])
    df["death_event"] = coerce_binary_series(df["death_event"])
    df["pfs_event"] = coerce_binary_series(df["pfs_event"])
    df["benefit"] = coerce_binary_series(df["benefit"], positive_values=benefit_cfg.get("positive_values"), negative_values=benefit_cfg.get("negative_values"))

    restricted_specs = derivation_cfg.get("restricted_survival", []) or []
    if isinstance(restricted_specs, dict):
        restricted_specs = [restricted_specs]
    for spec in restricted_specs:
        time_col = str(spec.get("time_col", ""))
        event_col = str(spec.get("event_col", ""))
        overwrite = bool(spec.get("overwrite", True))
        restricted_time, restricted_event = derive_restricted_survival(df, spec)
        if time_col not in df.columns or overwrite:
            df[time_col] = restricted_time
        else:
            df[time_col] = df[time_col].fillna(restricted_time)
        if event_col not in df.columns or overwrite:
            df[event_col] = restricted_event
        else:
            df[event_col] = coerce_binary_series(df[event_col]).fillna(restricted_event)
        df[event_col] = coerce_binary_series(df[event_col])

    binary_survival_specs = derivation_cfg.get("binary_survival", []) or []
    if isinstance(binary_survival_specs, dict):
        binary_survival_specs = [binary_survival_specs]
    for spec in binary_survival_specs:
        outcome_col = str(spec.get("outcome_col", ""))
        label_col = str(spec.get("label_col", ""))
        if not outcome_col:
            raise ValueError("Each clinical.derivation.binary_survival entry must define outcome_col.")
        overwrite = bool(spec.get("overwrite", True))
        binary_outcome, binary_labels = derive_binary_survival_endpoint(df, spec)
        if outcome_col not in df.columns or overwrite:
            df[outcome_col] = binary_outcome
        else:
            df[outcome_col] = coerce_binary_series(df[outcome_col]).fillna(binary_outcome)
        df[outcome_col] = coerce_binary_series(df[outcome_col])
        if label_col:
            if label_col not in df.columns or overwrite:
                df[label_col] = binary_labels
            else:
                df[label_col] = df[label_col].astype("string").fillna(binary_labels)

    return df.drop_duplicates(subset=["patient_id"]).reset_index(drop=True)


def load_clinical_df(cfg: dict) -> pd.DataFrame:
    clinical_cfg = cfg["clinical"]
    raw = read_table(clinical_cfg["path"], sheet_name=clinical_cfg.get("sheet_name", 0))
    if raw is None:
        raise FileNotFoundError(f"Clinical file not found: {clinical_cfg['path']}")
    return prepare_clinical_df(raw, clinical_cfg["column_map"], clinical_cfg.get("derivation", {}))


def read_selected_patients(path: str) -> set[str]:
    path_obj = Path(path)
    if not path_obj.exists():
        raise FileNotFoundError(f"Selected-patients file not found: {path}")

    if path_obj.suffix.lower() in {".txt", ".tsv", ".csv"}:
        lines = [normalize_id(line) for line in path_obj.read_text(encoding="utf-8").splitlines() if line.strip()]
        if lines and all("\t" not in line and "," not in line for line in lines):
            return set(lines)

    selected = read_table(path)
    if selected is None:
        raise FileNotFoundError(f"Selected-patients file not found: {path}")
    first_col = selected.columns[0]
    values = selected[first_col].dropna().tolist()
    values = [first_col, *values] if pd.notna(first_col) else values
    return {normalize_id(value) for value in values if str(value).strip()}


def infer_feature_family(feature_name: str) -> str:
    feature = str(feature_name)
    stem = ""
    core = feature
    if "__" in feature:
        left, right = feature.split("__", 1)
        known_prefixes = ("Centroid_", "curated_gene_", "denovo_gene_", "ilr_", "cellrank_", "cellrank_cr_b3")
        if not left.startswith(known_prefixes):
            stem = left
            core = right
    if core.startswith("Centroid_"):
        return f"centroid::{stem}" if stem else "centroid"
    if core.startswith("curated_gene_"):
        suffixes = ("__eigengene", "__pc1", "__pc2", "__singscore", "__ssgsea2_es")
        signature_stem = core
        for suffix in suffixes:
            if core.endswith(suffix):
                signature_stem = core[: -len(suffix)]
                break
        return f"signature::{stem}::{signature_stem}" if stem else f"signature::{signature_stem}"
    if core.startswith("cellrank_"):
        return f"cellrank::{stem}" if stem else "cellrank"
    if stem:
        return f"file::{stem}"
    if feature.startswith("PC"):
        return "pc"
    return "misc"


def infer_family_cap_group(feature_name: str) -> Optional[str]:
    feature = str(feature_name)
    stem = ""
    core = feature
    if "__" in feature:
        left, right = feature.split("__", 1)
        known_prefixes = ("Centroid_", "curated_gene_", "denovo_gene_", "ilr_", "cellrank_", "cellrank_cr_b3")
        if not left.startswith(known_prefixes):
            stem = left
            core = right
    if not core.startswith("curated_gene_"):
        return None

    suffixes = ("__eigengene", "__pc1", "__pc2", "__singscore", "__ssgsea2_es")
    signature_stem = core
    for suffix in suffixes:
        if core.endswith(suffix):
            signature_stem = core[: -len(suffix)]
            break
    return f"signature::{stem}::{signature_stem}" if stem else f"signature::{signature_stem}"


def infer_feature_level(feature_name: str) -> str:
    feature = str(feature_name)
    core = feature
    if "__" in feature:
        left, right = feature.split("__", 1)
        known_prefixes = ("Centroid_", "curated_gene_", "denovo_gene_", "ilr_", "cellrank_", "cellrank_cr_b3")
        core = feature if left.startswith(known_prefixes) else right
    if core.startswith("Centroid_"):
        return "PCA_centroid_continuum"
    if core.startswith("curated_gene_") and core.endswith("__ssgsea2_es"):
        return "pseudobulk_curated_ssgsea2"
    if core.startswith("curated_gene_") and core.endswith("__singscore"):
        return "pseudobulk_curated_singscore"
    if core.startswith("curated_gene_") and core.endswith("__pc1"):
        return "pseudobulk_curated_pc1"
    if core.startswith("curated_gene_") and core.endswith("__pc2"):
        return "pseudobulk_curated_pc2"
    if core.startswith("curated_gene_") and core.endswith("__eigengene"):
        return "pseudobulk_curated_eigengene"
    if core.startswith("denovo_gene_"):
        return "pseudobulk_de_novo"
    if core.startswith("ilr_"):
        return "ILR_composition"
    if core.startswith("cellrank_cr_b3"):
        return "status_dynamics"
    return "patient_feature"


def drop_bad_columns(df: pd.DataFrame, columns: list[str], min_non_missing_fraction: float, min_unique_values: int) -> list[str]:
    keep: list[str] = []
    n_rows = max(len(df), 1)
    for column in columns:
        non_missing_fraction = df[column].notna().sum() / n_rows
        n_unique = df[column].dropna().nunique()
        if non_missing_fraction >= min_non_missing_fraction and n_unique >= min_unique_values:
            keep.append(column)
    return keep


def _feature_input_specs(feature_cfg: dict) -> list[dict]:
    specs: list[dict] = []
    if feature_cfg.get("tables"):
        for item in feature_cfg.get("tables", []):
            if isinstance(item, str):
                specs.append({"path": item})
            else:
                specs.append(dict(item))
    elif feature_cfg.get("paths"):
        for path in feature_cfg.get("paths", []):
            specs.append({"path": path})
    elif feature_cfg.get("path"):
        specs.append({"path": feature_cfg.get("path")})
    else:
        raise ValueError("feature_library must provide one of: path, paths, or tables")
    return specs


def _make_prefix(spec: dict, path: str, idx: int, multiple: bool, feature_cfg: dict) -> str:
    mode = str(feature_cfg.get("prefix_mode", "auto"))
    if mode == "none":
        return ""
    alias = spec.get("name") or spec.get("alias")
    if alias:
        return str(alias)
    if mode == "auto" and not multiple:
        return ""
    return Path(path).stem


def _load_one_feature_table(spec: dict, feature_cfg: dict, selected: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    path = spec["path"]
    raw = read_table(path, sheet_name=spec.get("sheet_name", feature_cfg.get("sheet_name", 0)))
    if raw is None:
        raise FileNotFoundError(f"Feature library not found: {path}")

    id_col = str(spec.get("id_col", feature_cfg.get("id_col", "patient_id")))
    if id_col not in raw.columns:
        raise ValueError(f"Feature library is missing id column '{id_col}': {path}")

    raw = raw.copy()
    raw[id_col] = raw[id_col].map(normalize_id)
    if selected:
        raw = raw[raw[id_col].isin(selected)].copy()

    filters = spec.get("filters", feature_cfg.get("filters", {})) or {}
    for column, expected in filters.items():
        if column not in raw.columns:
            raise ValueError(f"Configured feature_library.filter column not found in {path}: {column}")
        allowed = expected if isinstance(expected, list) else [expected]
        raw = raw[raw[column].astype(str).isin({str(x) for x in allowed})].copy()

    carry_columns = [column for column in feature_cfg.get("carry_columns", []) if column in raw.columns]
    drop_columns = {id_col, *carry_columns, *feature_cfg.get("drop_columns", [])}

    include_regex = spec.get("candidate_include_regex", feature_cfg.get("candidate_include_regex"))
    exclude_regex = spec.get("candidate_exclude_regex", feature_cfg.get("candidate_exclude_regex"))
    include_re = re.compile(include_regex) if include_regex else None
    exclude_re = re.compile(exclude_regex) if exclude_regex else None

    candidate_features: list[str] = []
    for column in raw.columns:
        if column in drop_columns:
            continue
        if include_re and not include_re.search(str(column)):
            continue
        if exclude_re and exclude_re.search(str(column)):
            continue
        converted = pd.to_numeric(raw[column], errors="coerce")
        if converted.notna().sum() == 0:
            continue
        raw[column] = converted
        candidate_features.append(str(column))

    trace_df = raw[[id_col] + carry_columns].copy().rename(columns={id_col: "patient_id"})
    feature_df = raw[[id_col] + candidate_features].copy().rename(columns={id_col: "patient_id"})
    manifest = {
        "path": str(path),
        "n_rows": int(raw.shape[0]),
        "n_candidate_features": int(len(candidate_features)),
        "carry_columns": carry_columns,
    }
    return feature_df, trace_df, candidate_features, manifest


def _merge_trace_frames(trace_frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not trace_frames:
        return pd.DataFrame(columns=["patient_id"])
    merged = trace_frames[0].copy()
    for df in trace_frames[1:]:
        overlap = [c for c in df.columns if c in merged.columns and c != "patient_id"]
        if overlap:
            renamed = df.rename(columns={c: f"{c}__dup" for c in overlap})
            merged = merged.merge(renamed, on="patient_id", how="outer")
            for c in overlap:
                dup = f"{c}__dup"
                merged[c] = merged[c].combine_first(merged[dup])
                merged = merged.drop(columns=[dup])
        else:
            merged = merged.merge(df, on="patient_id", how="outer")
    return merged


def load_feature_library(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    feature_cfg = cfg["feature_library"]
    specs = _feature_input_specs(feature_cfg)
    selected = read_selected_patients(feature_cfg["selected_patients_file"]) if feature_cfg.get("selected_patients_file") else set()

    feature_frames: list[pd.DataFrame] = []
    trace_frames: list[pd.DataFrame] = []
    catalog_rows: list[dict] = []
    manifests: list[dict] = []

    multiple = len(specs) > 1
    for idx, spec in enumerate(specs, start=1):
        feature_df, trace_df, candidate_features, manifest = _load_one_feature_table(spec, feature_cfg, selected)
        prefix = _make_prefix(spec, spec["path"], idx, multiple, feature_cfg)

        rename_map = {}
        for feature in candidate_features:
            new_name = f"{prefix}__{feature}" if prefix else feature
            rename_map[feature] = new_name
            catalog_rows.append({
                "feature_name": new_name,
                "source_feature_name": feature,
                "source_path": str(spec["path"]),
                "source_alias": prefix or Path(spec["path"]).stem,
            })
        feature_df = feature_df.rename(columns=rename_map)
        feature_frames.append(feature_df)
        trace_frames.append(trace_df)
        manifest["prefix"] = prefix
        manifests.append(manifest)

    merged_feature_df: Optional[pd.DataFrame] = None
    for df in feature_frames:
        merged_feature_df = df if merged_feature_df is None else merged_feature_df.merge(df, on="patient_id", how="outer")
    assert merged_feature_df is not None

    merged_trace_df = _merge_trace_frames(trace_frames)
    candidate_features = [c for c in merged_feature_df.columns if c != "patient_id"]
    if not candidate_features:
        raise ValueError("No usable numeric candidate features were found across feature tables.")

    catalog = pd.DataFrame(catalog_rows)
    if feature_cfg.get("feature_catalog_path"):
        extra_catalog = read_table(feature_cfg["feature_catalog_path"])
        if extra_catalog is None:
            raise FileNotFoundError(f"Feature catalog not found: {feature_cfg['feature_catalog_path']}")
        feature_col = str(feature_cfg.get("feature_catalog_feature_col", "feature_name"))
        if feature_col not in extra_catalog.columns:
            raise ValueError(f"Feature catalog is missing feature name column: {feature_col}")
        extra_catalog = extra_catalog.rename(columns={feature_col: "feature_name"})
        catalog = catalog.merge(extra_catalog, on="feature_name", how="left")

    if "family" not in catalog.columns:
        catalog["family"] = catalog["feature_name"].map(infer_feature_family)
    else:
        catalog["family"] = catalog["family"].fillna(catalog["feature_name"].map(infer_feature_family))
    if "family_cap_group" not in catalog.columns:
        catalog["family_cap_group"] = catalog["feature_name"].map(infer_family_cap_group)
    else:
        catalog["family_cap_group"] = catalog["family_cap_group"].fillna(catalog["feature_name"].map(infer_family_cap_group))
    catalog["family_cap_applies"] = catalog["family_cap_group"].notna()
    if "feature_level" not in catalog.columns:
        catalog["feature_level"] = catalog["feature_name"].map(infer_feature_level)
    else:
        catalog["feature_level"] = catalog["feature_level"].fillna(catalog["feature_name"].map(infer_feature_level))
    if "description" not in catalog.columns:
        catalog["description"] = catalog["feature_name"].map(lambda x: str(x).replace("_", " "))
    else:
        catalog["description"] = catalog["description"].fillna(catalog["feature_name"].map(lambda x: str(x).replace("_", " ")))

    catalog = catalog.drop_duplicates(subset=["feature_name"]).reset_index(drop=True)
    manifest = {
        "feature_library_inputs": manifests,
        "n_feature_tables": len(specs),
        "n_feature_rows": int(merged_feature_df.shape[0]),
        "n_candidate_features": int(len(candidate_features)),
        "selected_patients_file": feature_cfg.get("selected_patients_file"),
        "n_selected_patients_filter": int(len(selected)),
        "carry_columns": [c for c in merged_trace_df.columns if c != "patient_id"],
    }
    return merged_feature_df, merged_trace_df, manifest, catalog


def build_search_dataset(cfg: dict) -> SearchDataset:
    clinical_df = load_clinical_df(cfg)
    feature_df, trace_df, feature_manifest, feature_catalog = load_feature_library(cfg)

    merged = clinical_df.merge(trace_df, on="patient_id", how="inner")
    merged = merged.merge(feature_df, on="patient_id", how="inner")
    if merged.empty:
        raise ValueError("No matched patients remain after merging clinical data and feature library.")

    carry_columns = [column for column in trace_df.columns if column != "patient_id" and column in merged.columns]
    clinical_feature_sets = resolve_clinical_feature_sets(cfg, merged.columns.tolist())
    raw_candidate_features = [feature for feature in feature_df.columns if feature != "patient_id" and feature in merged.columns]

    qc_cfg = cfg.get("modeling", {}).get("data_qc", {})
    min_non_missing_fraction = float(qc_cfg.get("min_non_missing_fraction", 0.6))
    min_unique_values = int(qc_cfg.get("min_unique_values", 2))
    candidate_features = drop_bad_columns(
        merged,
        columns=raw_candidate_features,
        min_non_missing_fraction=min_non_missing_fraction,
        min_unique_values=min_unique_values,
    )
    if not candidate_features:
        raise ValueError("No candidate immune features remain after applying data_qc filters.")

    feature_catalog = feature_catalog.copy()
    feature_catalog["non_missing_fraction"] = feature_catalog["feature_name"].map(
        lambda col: float(merged[col].notna().mean()) if col in merged.columns else np.nan
    )
    feature_catalog["n_unique"] = feature_catalog["feature_name"].map(
        lambda col: int(merged[col].dropna().nunique()) if col in merged.columns else 0
    )
    feature_catalog["kept_for_modeling"] = feature_catalog["feature_name"].isin(candidate_features)
    feature_catalog = feature_catalog.sort_values(
        ["kept_for_modeling", "feature_level", "feature_name"],
        ascending=[False, True, True],
    ).reset_index(drop=True)

    manifest = {
        "clinical_path": str(cfg["clinical"]["path"]),
        "n_clinical_rows": int(clinical_df.shape[0]),
        "n_matched_patients": int(merged.shape[0]),
        "data_qc": {
            "min_non_missing_fraction": min_non_missing_fraction,
            "min_unique_values": min_unique_values,
        },
        "n_candidate_features_before_qc": int(len(raw_candidate_features)),
        "n_candidate_features_after_qc": int(len(candidate_features)),
        **feature_manifest,
    }

    return SearchDataset(
        analysis_df=merged,
        clinical_feature_sets=clinical_feature_sets,
        candidate_features=candidate_features,
        carry_columns=carry_columns,
        feature_catalog=feature_catalog,
        manifest=manifest,
    )
