"""Prepare a matched patient-level modeling table from clinical and immune inputs."""

from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd
from scipy.linalg import helmert

from .design import FeatureBlock, PreparedDataset, resolve_clinical_feature_sets
from .io_utils import make_unique_names, normalize_id, read_table, sanitize_feature_name


ESSENTIAL_CLINICAL_COLUMNS = {
    "patient_id",
    "sex",
    "age",
    "histology",
    "smoking",
    "ecog",
    "drug",
    "pd_event",
    "pfs_days",
    "egfr_status",
    "io_line",
    "previous_chemo",
    "previous_target",
    "pd_l1_tps",
}
OPTIONAL_CLINICAL_COLUMNS = {"death_event", "pfs_event", "recist", "benefit"}


def infer_composition_representation(values: np.ndarray) -> str:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return "raw"
    row_sums = np.nansum(values, axis=1)
    frac_zero_sum = float(np.mean(np.isclose(row_sums, 0.0, atol=1e-6))) if row_sums.size else 0.0
    frac_negative = float(np.mean(finite < -1e-8))
    if frac_negative > 0.05 and frac_zero_sum > 0.5:
        return "clr"
    return "raw"


def close_composition(values: np.ndarray) -> np.ndarray:
    closed = np.asarray(values, dtype=float).copy()
    closed = np.clip(closed, a_min=0.0, a_max=None)
    row_sums = closed.sum(axis=1, keepdims=True)
    valid = row_sums.ravel() > 0
    if np.any(valid):
        closed[valid] = closed[valid] / row_sums[valid]
    if np.any(~valid):
        closed[~valid] = np.nan
    return closed


def clr_to_composition(values: np.ndarray) -> np.ndarray:
    exp_values = np.exp(np.asarray(values, dtype=float))
    row_sums = exp_values.sum(axis=1, keepdims=True)
    valid = row_sums.ravel() > 0
    out = np.full_like(exp_values, np.nan, dtype=float)
    if np.any(valid):
        out[valid] = exp_values[valid] / row_sums[valid]
    return out


def composition_to_clr(values: np.ndarray) -> np.ndarray:
    log_prop = np.log(np.asarray(values, dtype=float))
    return log_prop - log_prop.mean(axis=1, keepdims=True)


def composition_to_ilr(values: np.ndarray) -> np.ndarray:
    n_parts = values.shape[1]
    if n_parts < 2:
        return np.empty((values.shape[0], 0), dtype=float)
    basis = helmert(n_parts, full=False).T
    clr_values = composition_to_clr(values)
    return clr_values @ basis


def resolve_compositional_input(
    values: np.ndarray,
    input_representation: str,
) -> tuple[np.ndarray, str]:
    representation = str(input_representation or "auto").strip().lower()
    if representation == "auto":
        representation = infer_composition_representation(values)

    if representation in {"raw", "closed", "proportion", "proportions"}:
        return close_composition(values), "raw"
    if representation == "clr":
        return clr_to_composition(values), "clr"
    raise ValueError(
        "Unsupported compositional input representation. "
        "Use one of: auto, raw, closed, proportion, proportions, clr."
    )


def extract_pc_index(column_name: str) -> Optional[int]:
    match = re.search(r"_PC(\d+)$", str(column_name))
    if not match:
        return None
    return int(match.group(1))


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
        "pembrolizumab",
        "nivolumab",
        "cemiplimab",
        "dostarlimab",
        "tislelizumab",
        "toripalimab",
        "sintilimab",
        "camrelizumab",
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


def coerce_binary_series(
    series: pd.Series,
    positive_values: Optional[list[object]] = None,
    negative_values: Optional[list[object]] = None,
) -> pd.Series:
    out = pd.Series(np.nan, index=series.index, dtype=float)
    numeric = pd.to_numeric(series, errors="coerce")
    mask_numeric = numeric.notna()
    if mask_numeric.any():
        out.loc[mask_numeric & (numeric == 1)] = 1.0
        out.loc[mask_numeric & (numeric == 0)] = 0.0

    positive_set = {
        str(value).strip().upper()
        for value in (positive_values or [1, "1", "YES", "Y", "TRUE", "BENEFIT", "RESPONDER", "DCB"])
    }
    negative_set = {
        str(value).strip().upper()
        for value in (negative_values or [0, "0", "NO", "N", "FALSE", "NO BENEFIT", "NONBENEFIT", "NON-BENEFIT", "NCB"])
    }
    text = series.astype("string").str.strip().str.upper()

    out.loc[text.isin(positive_set)] = 1.0
    out.loc[text.isin(negative_set)] = 0.0
    return out


def derive_pfs_event(df: pd.DataFrame, derivation_cfg: dict) -> pd.Series:
    pfs_cfg = derivation_cfg.get("pfs_event", {})
    method = str(pfs_cfg.get("method", "any_of")).lower()
    source_columns = [str(column) for column in pfs_cfg.get("source_columns", ["pd_event", "death_event"])]
    usable = [column for column in source_columns if column in df.columns]
    if not usable:
        return pd.Series(np.nan, index=df.index, dtype=float)

    source_df = pd.concat([coerce_binary_series(df[column]) for column in usable], axis=1)
    if method == "any_of":
        return source_df.max(axis=1, skipna=True)
    if method == "all_of":
        return source_df.min(axis=1, skipna=True)
    raise ValueError(f"Unsupported clinical.derivation.pfs_event.method: {method}")


def derive_benefit(df: pd.DataFrame, derivation_cfg: dict) -> pd.Series:
    benefit_cfg = derivation_cfg.get("benefit", {})
    recist_column = str(benefit_cfg.get("recist_column", "recist"))
    time_column = str(benefit_cfg.get("time_column", "pfs_days"))
    if recist_column not in df.columns or time_column not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)

    positive_response_values = {
        str(value).strip().upper() for value in benefit_cfg.get("positive_response_values", ["CR", "PR"])
    }
    durable_stable_values = {
        str(value).strip().upper() for value in benefit_cfg.get("durable_stable_values", ["SD"])
    }
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
    # Restrict follow-up to a fixed horizon while preserving time-to-event information before that horizon.
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


def prepare_clinical_df(raw: pd.DataFrame, column_map: dict[str, Optional[str]], derivation_cfg: Optional[dict] = None) -> pd.DataFrame:
    # Standardize raw clinical metadata columns into a stable modeling schema.
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
    benefit_existing = coerce_binary_series(
        df["benefit"],
        positive_values=benefit_cfg.get("positive_values"),
        negative_values=benefit_cfg.get("negative_values"),
    )
    prefer_existing = bool(benefit_cfg.get("prefer_existing", True))
    derived_benefit = derive_benefit(df, derivation_cfg)
    if prefer_existing:
        df["benefit"] = benefit_existing.fillna(derived_benefit)
    else:
        df["benefit"] = derived_benefit.fillna(benefit_existing)

    df["pd_event"] = coerce_binary_series(df["pd_event"])
    df["death_event"] = coerce_binary_series(df["death_event"])
    df["pfs_event"] = coerce_binary_series(df["pfs_event"])
    df["benefit"] = coerce_binary_series(
        df["benefit"],
        positive_values=benefit_cfg.get("positive_values"),
        negative_values=benefit_cfg.get("negative_values"),
    )

    restricted_specs = derivation_cfg.get("restricted_survival", []) or []
    if isinstance(restricted_specs, dict):
        restricted_specs = [restricted_specs]
    for spec in restricted_specs:
        time_col = str(spec.get("time_col", ""))
        event_col = str(spec.get("event_col", ""))
        if not time_col or not event_col:
            raise ValueError("Each clinical.derivation.restricted_survival entry must define time_col and event_col.")
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

    df = df.drop_duplicates(subset=["patient_id"]).reset_index(drop=True)
    return df


def load_clinical_df(cfg: dict) -> pd.DataFrame:
    clinical_cfg = cfg["clinical"]
    raw = read_table(clinical_cfg["path"], sheet_name=clinical_cfg.get("sheet_name", 0))
    if raw is None:
        raise FileNotFoundError(f"Clinical table not found: {clinical_cfg['path']}")
    return prepare_clinical_df(raw, clinical_cfg["column_map"], clinical_cfg.get("derivation", {}))


def apply_row_filters(df: pd.DataFrame, filters: dict[str, object]) -> pd.DataFrame:
    out = df.copy()
    for column, expected in filters.items():
        if column not in out.columns:
            raise ValueError(f"Configured immune_features.filter column not found: {column}")
        allowed = expected if isinstance(expected, list) else [expected]
        allowed_text = {str(value) for value in allowed}
        out = out[out[column].astype(str).isin(allowed_text)].copy()
    return out


def resolve_feature_table_path(immune_cfg: dict) -> str:
    path = immune_cfg.get("feature_table_path") or immune_cfg.get("path")
    if not path:
        raise ValueError("immune_features.feature_table_path or immune_features.path must be provided.")
    return str(path)


def build_proportion_block_from_long(
    long_df: pd.DataFrame,
    trace_df: pd.DataFrame,
    block_name: str,
    block_cfg: dict,
    trace_analysis_id_col: str,
) -> tuple[Optional[pd.DataFrame], list[dict[str, object]]]:
    # Reconstruct raw composition features directly from long-form subtype outputs, then transform for modeling.
    analysis_id_col = str(block_cfg.get("analysis_id_col", "analysis_id"))
    celltype_col = str(block_cfg.get("celltype_col", "celltype"))
    value_col = str(block_cfg.get("value_col", "prop"))
    count_col = block_cfg.get("count_col")
    prefix = str(block_cfg.get("prefix", block_name))
    transform = str(block_cfg.get("transform", "none")).lower()
    required = bool(block_cfg.get("required", False))
    include_celltypes = [str(value) for value in block_cfg.get("include_celltypes", [])]

    required_columns = {analysis_id_col, celltype_col, value_col}
    missing_columns = [column for column in required_columns if column not in long_df.columns]
    if missing_columns:
        raise ValueError(
            f"Proportion-long source for block '{block_name}' is missing required columns: {missing_columns}"
        )
    if trace_analysis_id_col not in trace_df.columns:
        raise ValueError(
            f"Trace table is missing analysis-id column '{trace_analysis_id_col}' needed for block '{block_name}'."
        )

    long_use = long_df.copy()
    long_use[analysis_id_col] = long_use[analysis_id_col].astype(str)
    trace_use = trace_df[["patient_id", trace_analysis_id_col]].copy()
    trace_use[trace_analysis_id_col] = trace_use[trace_analysis_id_col].astype(str)
    long_use = long_use[long_use[analysis_id_col].isin(trace_use[trace_analysis_id_col])].copy()

    if include_celltypes:
        long_use = long_use[long_use[celltype_col].astype(str).isin(include_celltypes)].copy()
        celltype_order = include_celltypes
    else:
        celltype_order = sorted(long_use[celltype_col].dropna().astype(str).unique().tolist())

    if long_use.empty or not celltype_order:
        if required:
            raise ValueError(f"Required immune block '{block_name}' has no usable rows in the proportion-long source.")
        return None, []

    if count_col and str(count_col) in long_use.columns:
        count_column = str(count_col)
        count_wide = (
            long_use.pivot_table(index=analysis_id_col, columns=celltype_col, values=count_column, aggfunc="sum")
            .reindex(index=trace_use[trace_analysis_id_col], columns=celltype_order)
            .fillna(0.0)
        )
        pseudocount = float(block_cfg.get("pseudocount", 0.0))
        if pseudocount < 0:
            raise ValueError(f"Immune block '{block_name}' requires pseudocount >= 0.")
        composition_df = (count_wide + pseudocount).div((count_wide + pseudocount).sum(axis=1), axis=0)
        source_representation = "count_smoothed_composition"
    else:
        value_wide = (
            long_use.pivot_table(index=analysis_id_col, columns=celltype_col, values=value_col, aggfunc="sum")
            .reindex(index=trace_use[trace_analysis_id_col], columns=celltype_order)
            .fillna(0.0)
        )
        composition_df = pd.DataFrame(
            close_composition(value_wide.to_numpy(dtype=float)),
            index=value_wide.index,
            columns=value_wide.columns,
        )
        source_representation = "raw_composition"

    block = pd.DataFrame({"patient_id": trace_use["patient_id"].to_numpy()})
    feature_rows: list[dict[str, object]] = []
    if transform == "ilr":
        transformed = composition_to_ilr(composition_df.to_numpy(dtype=float))
        if transformed.shape[1] == 0:
            if required:
                raise ValueError(f"Immune block '{block_name}' needs at least two cell types for ILR transform.")
            return None, []
        unique_names = make_unique_names([f"{prefix}__ilr{i + 1}" for i in range(transformed.shape[1])])
        transformed_df = pd.DataFrame(transformed, index=composition_df.index, columns=unique_names)
        for column in unique_names:
            block[column] = transformed_df[column].to_numpy()
            feature_rows.append(
                {
                    "block": block_name,
                    "source_column": "__ilr_balance__",
                    "source_columns": "|".join(celltype_order),
                    "feature_column": column,
                    "transform": transform,
                    "input_representation": source_representation,
                    "prefix": prefix,
                }
            )
    elif transform == "clr":
        transformed_df = pd.DataFrame(
            composition_to_clr(composition_df.to_numpy(dtype=float)),
            index=composition_df.index,
            columns=celltype_order,
        )
        unique_names = make_unique_names([f"{prefix}__{sanitize_feature_name(column)}" for column in celltype_order])
        for source_column, feature_column in zip(celltype_order, unique_names):
            block[feature_column] = transformed_df[source_column].to_numpy()
            feature_rows.append(
                {
                    "block": block_name,
                    "source_column": source_column,
                    "source_columns": source_column,
                    "feature_column": feature_column,
                    "transform": transform,
                    "input_representation": source_representation,
                    "prefix": prefix,
                }
            )
    elif transform in {"none", "", "raw"}:
        unique_names = make_unique_names([f"{prefix}__{sanitize_feature_name(column)}" for column in celltype_order])
        for source_column, feature_column in zip(celltype_order, unique_names):
            block[feature_column] = composition_df[source_column].to_numpy()
            feature_rows.append(
                {
                    "block": block_name,
                    "source_column": source_column,
                    "source_columns": source_column,
                    "feature_column": feature_column,
                    "transform": "raw",
                    "input_representation": source_representation,
                    "prefix": prefix,
                }
            )
    else:
        raise ValueError(
            f"Unsupported immune transform '{transform}' for block '{block_name}' from proportion-long source."
        )
    return block, feature_rows


def extract_immune_block(
    immune_df: pd.DataFrame,
    patient_id_col: str,
    block_name: str,
    block_cfg: dict,
    protected_columns: set[str],
) -> tuple[Optional[pd.DataFrame], list[dict[str, object]]]:
    # Build one immune feature block at a time so proportions, PCs, and future blocks stay modular.
    include_columns = [col for col in block_cfg.get("include_columns", []) if col in immune_df.columns]
    exclude_columns = set(block_cfg.get("exclude_columns", []))
    include_regex = block_cfg.get("include_regex")
    exclude_regex = block_cfg.get("exclude_regex")
    prefix = str(block_cfg.get("prefix", block_name))
    transform = str(block_cfg.get("transform", "none")).lower()
    input_representation = str(block_cfg.get("input_representation", "auto"))
    required = bool(block_cfg.get("required", False))

    if include_columns:
        candidates = include_columns
    else:
        candidates = [col for col in immune_df.columns if col not in protected_columns]

    if include_regex:
        regex = re.compile(include_regex)
        candidates = [col for col in candidates if regex.search(col)]
    if exclude_regex:
        regex = re.compile(exclude_regex)
        candidates = [col for col in candidates if not regex.search(col)]

    include_pc_indices = {int(value) for value in block_cfg.get("include_pc_indices", [])}
    max_pc_index = block_cfg.get("max_pc_index")
    max_pc_index = int(max_pc_index) if max_pc_index is not None else None
    if include_pc_indices or max_pc_index is not None:
        filtered_candidates = []
        for column in candidates:
            pc_index = extract_pc_index(column)
            if pc_index is None:
                continue
            if include_pc_indices and pc_index not in include_pc_indices:
                continue
            if max_pc_index is not None and pc_index > max_pc_index:
                continue
            filtered_candidates.append(column)
        candidates = filtered_candidates

    candidates = [col for col in candidates if col not in exclude_columns and col != patient_id_col]
    numeric_columns = []
    converted_cache: dict[str, pd.Series] = {}
    for column in candidates:
        converted = pd.to_numeric(immune_df[column], errors="coerce")
        if converted.notna().sum() > 0:
            converted_cache[column] = converted
            numeric_columns.append(column)

    if not numeric_columns:
        if required:
            raise ValueError(f"Required immune block '{block_name}' has no usable numeric columns.")
        return None, []

    source_matrix = pd.DataFrame({column: converted_cache[column] for column in numeric_columns})
    block = pd.DataFrame({"patient_id": immune_df[patient_id_col]})
    feature_rows: list[dict[str, object]] = []

    if transform in {"clr", "ilr"}:
        pseudocount = float(block_cfg.get("pseudocount", 1e-6))
        if pseudocount <= 0:
            raise ValueError(f"Immune block '{block_name}' requires pseudocount > 0 for {transform} transform.")
        source_values = source_matrix.to_numpy(dtype=float)
        composition, resolved_input = resolve_compositional_input(source_values, input_representation)
        adjusted = close_composition(np.nan_to_num(composition, nan=0.0) + pseudocount)
        if transform == "clr":
            transformed_values = composition_to_clr(adjusted)
            generated_names = [f"{prefix}__{sanitize_feature_name(column)}" for column in numeric_columns]
            unique_names = make_unique_names(generated_names)
            for source_column, feature_column in zip(numeric_columns, unique_names):
                feature_rows.append(
                    {
                        "block": block_name,
                        "source_column": source_column,
                        "source_columns": source_column,
                        "feature_column": feature_column,
                        "transform": transform,
                        "input_representation": resolved_input,
                        "prefix": prefix,
                    }
                )
        else:
            transformed_values = composition_to_ilr(adjusted)
            if transformed_values.shape[1] == 0:
                if required:
                    raise ValueError(f"Immune block '{block_name}' needs at least two source columns for ILR transform.")
                return None, []
            generated_names = [f"{prefix}__ilr{i + 1}" for i in range(transformed_values.shape[1])]
            unique_names = make_unique_names(generated_names)
            source_columns_joined = "|".join(numeric_columns)
            for feature_column in unique_names:
                feature_rows.append(
                    {
                        "block": block_name,
                        "source_column": "__ilr_balance__",
                        "source_columns": source_columns_joined,
                        "feature_column": feature_column,
                        "transform": transform,
                        "input_representation": resolved_input,
                        "prefix": prefix,
                    }
                )
        transformed_df = pd.DataFrame(transformed_values, index=source_matrix.index, columns=unique_names)
        for column in unique_names:
            block[column] = transformed_df[column]
    elif transform in {"none", "", "raw"}:
        generated_names = [f"{prefix}__{sanitize_feature_name(column)}" for column in numeric_columns]
        unique_names = make_unique_names(generated_names)
        rename_map = dict(zip(numeric_columns, unique_names))
        for source_column, feature_column in zip(numeric_columns, unique_names):
            block[feature_column] = source_matrix[source_column]
            feature_rows.append(
                {
                    "block": block_name,
                    "source_column": source_column,
                    "source_columns": source_column,
                    "feature_column": feature_column,
                    "transform": "none",
                    "input_representation": input_representation,
                    "prefix": prefix,
                }
            )
    else:
        raise ValueError(f"Unsupported immune transform '{transform}' for block '{block_name}'.")

    return block, feature_rows


def load_immune_blocks(cfg: dict) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame, dict[str, pd.DataFrame], dict[str, object]]:
    # Read the feature table for patient trace/PC blocks and optionally separate long-form proportion tables.
    immune_cfg = cfg["immune_features"]
    feature_table_path = resolve_feature_table_path(immune_cfg)
    raw = read_table(feature_table_path, sheet_name=immune_cfg.get("sheet_name", 0))
    if raw is None:
        raise FileNotFoundError(f"Immune feature table not found: {feature_table_path}")

    if immune_cfg.get("filters"):
        raw = apply_row_filters(raw, immune_cfg["filters"])

    patient_id_col = immune_cfg.get("id_col", "patient_id")
    if patient_id_col not in raw.columns:
        raise ValueError(f"Immune feature table is missing id column: {patient_id_col}")

    raw = raw.copy()
    raw[patient_id_col] = raw[patient_id_col].map(normalize_id)
    if raw[patient_id_col].duplicated().any():
        dup_ids = raw.loc[raw[patient_id_col].duplicated(), patient_id_col].astype(str).tolist()[:10]
        raise ValueError(
            "Immune feature table must have one row per patient after filtering. "
            f"Duplicated ids include: {dup_ids}"
        )

    carry_columns = [column for column in immune_cfg.get("carry_columns", []) if column in raw.columns]
    trace_df = raw[[patient_id_col] + carry_columns].copy().rename(columns={patient_id_col: "patient_id"})
    diagnostic_tables: dict[str, pd.DataFrame] = {}
    analysis_id_col = str(immune_cfg.get("analysis_id_col", "analysis_id"))

    protected_columns = set(carry_columns) | set(immune_cfg.get("drop_columns", [])) | {patient_id_col}
    block_tables: dict[str, pd.DataFrame] = {}
    feature_rows: list[dict[str, object]] = []
    block_sources: dict[str, str] = {}
    for block_name, block_cfg in immune_cfg.get("blocks", {}).items():
        source_kind = str(block_cfg.get("source", "feature_table")).strip().lower()
        if source_kind == "proportion_long":
            block_path = block_cfg.get("path")
            if not block_path:
                raise ValueError(f"Immune block '{block_name}' with source 'proportion_long' must define a path.")
            block_raw = read_table(str(block_path), sheet_name=block_cfg.get("sheet_name", 0))
            if block_raw is None:
                raise FileNotFoundError(f"Immune block source not found: {block_path}")
            block_df, rows = build_proportion_block_from_long(
                long_df=block_raw,
                trace_df=trace_df,
                block_name=block_name,
                block_cfg=block_cfg,
                trace_analysis_id_col=analysis_id_col,
            )
            block_sources[block_name] = str(block_path)
        elif source_kind == "feature_table":
            block_df, rows = extract_immune_block(
                immune_df=raw,
                patient_id_col=patient_id_col,
                block_name=block_name,
                block_cfg=block_cfg,
                protected_columns=protected_columns,
            )
            block_sources[block_name] = feature_table_path
        else:
            raise ValueError(
                f"Unsupported immune block source '{source_kind}' for block '{block_name}'. "
                "Use 'feature_table' or 'proportion_long'."
            )
        if block_df is not None:
            block_tables[block_name] = block_df
            feature_rows.extend(rows)

    if not block_tables:
        raise ValueError("No usable immune feature blocks were extracted from the immune feature table.")

    feature_dictionary = pd.DataFrame(feature_rows)
    manifest = {
        "immune_feature_table_path": feature_table_path,
        "immune_block_sources": block_sources,
    }
    return trace_df, block_tables, feature_dictionary, diagnostic_tables, manifest


def drop_bad_columns(df: pd.DataFrame, columns: list[str], min_non_missing_fraction: float, min_unique_values: int) -> list[str]:
    keep: list[str] = []
    n_rows = max(len(df), 1)
    for column in columns:
        non_missing_fraction = df[column].notna().sum() / n_rows
        n_unique = df[column].dropna().nunique()
        if non_missing_fraction >= min_non_missing_fraction and n_unique >= min_unique_values:
            keep.append(column)
    return keep


def extract_source_columns_from_dictionary(feature_dictionary: pd.DataFrame) -> list[str]:
    source_columns: list[str] = []
    for value in feature_dictionary.get("source_columns", pd.Series(dtype=str)).dropna().astype(str):
        for source_column in value.split("|"):
            source_column = source_column.strip()
            if source_column:
                source_columns.append(source_column)
    return list(dict.fromkeys(source_columns))


def build_prepared_dataset(cfg: dict) -> PreparedDataset:
    # Step 1: create one matched cohort shared by all model variants for fair comparison.
    clinical_df = load_clinical_df(cfg)
    trace_df, block_tables, feature_dictionary, diagnostic_tables, immune_manifest = load_immune_blocks(cfg)

    merged = clinical_df.merge(trace_df, on="patient_id", how="inner")
    for block_df in block_tables.values():
        merged = merged.merge(block_df, on="patient_id", how="left")

    merged = merged.dropna(subset=["patient_id"]).reset_index(drop=True)
    if merged.empty:
        raise ValueError("No matched patients remain after merging clinical and immune inputs.")

    qc_cfg = cfg.get("modeling", {}).get("data_qc", {})
    min_non_missing_fraction = float(qc_cfg.get("min_non_missing_fraction", 0.6))
    min_unique_values = int(qc_cfg.get("min_unique_values", 2))

    immune_blocks: dict[str, FeatureBlock] = {}
    keep_columns = ["patient_id"]
    carry_columns = [column for column in trace_df.columns if column != "patient_id" and column in merged.columns]
    keep_columns.extend(carry_columns)

    feature_dictionary = feature_dictionary.copy()
    feature_dictionary["non_missing_fraction"] = np.nan
    feature_dictionary["n_unique"] = np.nan
    feature_dictionary["kept_for_modeling"] = False

    for block_name, block_df in block_tables.items():
        # Step 2: drop sparse or constant immune features before model fitting.
        raw_columns = [column for column in block_df.columns if column != "patient_id"]
        keep_block_columns = drop_bad_columns(
            merged,
            columns=raw_columns,
            min_non_missing_fraction=min_non_missing_fraction,
            min_unique_values=min_unique_values,
        )
        block_feature_dict = feature_dictionary["block"] == block_name
        for feature_column in raw_columns:
            if feature_column not in merged.columns:
                continue
            mask = feature_dictionary["feature_column"] == feature_column
            feature_dictionary.loc[mask, "non_missing_fraction"] = merged[feature_column].notna().mean()
            feature_dictionary.loc[mask, "n_unique"] = merged[feature_column].dropna().nunique()
            feature_dictionary.loc[mask, "kept_for_modeling"] = feature_column in keep_block_columns

        kept_feature_dict = feature_dictionary.loc[block_feature_dict & feature_dictionary["kept_for_modeling"]]
        source_columns = extract_source_columns_from_dictionary(kept_feature_dict)
        immune_blocks[block_name] = FeatureBlock(
            name=block_name,
            columns=keep_block_columns,
            source_columns=source_columns,
            transform=str(cfg["immune_features"]["blocks"][block_name].get("transform", "none")).lower(),
            prefix=str(cfg["immune_features"]["blocks"][block_name].get("prefix", block_name)),
        )
        keep_columns.extend(keep_block_columns)

    keep_columns.extend([column for column in clinical_df.columns if column not in {"patient_id", "drug"}])
    keep_columns = [column for column in keep_columns if column in merged.columns]
    keep_columns = list(dict.fromkeys(keep_columns))
    analysis_df = merged[keep_columns].copy()

    # Step 3: register the named clinical covariate sets that analyses can request later.
    clinical_feature_sets = resolve_clinical_feature_sets(cfg, available_columns=analysis_df.columns.tolist())

    manifest = {
        "clinical_path": str(cfg["clinical"]["path"]),
        "clinical_exclude_from_modeling": list(cfg.get("clinical", {}).get("exclude_from_modeling", []) or []),
        "n_clinical_rows": int(clinical_df.shape[0]),
        "n_immune_rows": int(trace_df.shape[0]),
        "n_matched_patients": int(analysis_df.shape[0]),
        "n_pd_event_nonmissing": int(analysis_df["pd_event"].notna().sum()) if "pd_event" in analysis_df.columns else 0,
        "n_pd_events": int(analysis_df["pd_event"].fillna(0).sum()) if "pd_event" in analysis_df.columns else 0,
        "n_death_event_nonmissing": int(analysis_df["death_event"].notna().sum()) if "death_event" in analysis_df.columns else 0,
        "n_death_events": int(analysis_df["death_event"].fillna(0).sum()) if "death_event" in analysis_df.columns else 0,
        "n_pfs_event_nonmissing": int(analysis_df["pfs_event"].notna().sum()) if "pfs_event" in analysis_df.columns else 0,
        "n_pfs_events": int(analysis_df["pfs_event"].fillna(0).sum()) if "pfs_event" in analysis_df.columns else 0,
        "n_benefit_nonmissing": int(analysis_df["benefit"].notna().sum()) if "benefit" in analysis_df.columns else 0,
        "n_benefit_events": int(analysis_df["benefit"].fillna(0).sum()) if "benefit" in analysis_df.columns else 0,
        "immune_blocks": {name: len(block.columns) for name, block in immune_blocks.items()},
        "clinical_feature_sets": {name: len(columns) for name, columns in clinical_feature_sets.items()},
        **immune_manifest,
    }

    return PreparedDataset(
        analysis_df=analysis_df,
        clinical_feature_sets=clinical_feature_sets,
        immune_blocks=immune_blocks,
        feature_dictionary=feature_dictionary.sort_values(["block", "feature_column"]).reset_index(drop=True),
        carry_columns=carry_columns,
        diagnostic_tables=diagnostic_tables,
        manifest=manifest,
    )
