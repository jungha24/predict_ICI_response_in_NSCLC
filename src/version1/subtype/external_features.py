from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


SUPPORTED_TEXT_SUFFIXES = {".txt", ".tsv", ".csv"}
SUPPORTED_PARQUET_SUFFIXES = {".parquet", ".pq"}


def first_non_null(series: pd.Series):
    vals = series.dropna()
    if len(vals) == 0:
        return np.nan
    return vals.iloc[0]


def normalize_patient_id(value: object) -> object:
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    return text or np.nan


def read_external_feature_table(path: str) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix in SUPPORTED_PARQUET_SUFFIXES:
        return pd.read_parquet(path)
    if suffix == '.csv':
        return pd.read_csv(path)
    if suffix in {'.txt', '.tsv'}:
        return pd.read_csv(path, sep='	')
    return pd.read_csv(path, sep=None, engine='python')


def load_external_patient_features(patient_meta: pd.DataFrame, cfg: dict) -> Tuple[pd.DataFrame, Dict[str, object]]:
    # Optionally merge precomputed patient-level immune features into the subtype feature table.
    external_cfg = cfg.get('external_patient_features', {})
    base = patient_meta[['analysis_id', 'patient_id']].copy()
    manifest: Dict[str, object] = {
        'enabled': False,
        'path': None,
        'feature_columns': [],
    }
    if not bool(external_cfg.get('enabled', False)):
        return base[['analysis_id']].copy(), manifest

    path = str(external_cfg.get('path', '')).strip()
    if not path:
        raise ValueError('external_patient_features.path must be provided when enabled=true')

    raw = read_external_feature_table(path)
    id_col = str(external_cfg.get('id_col', 'patient_id'))
    if id_col not in raw.columns:
        raise ValueError(f'external_patient_features.id_col not found: {id_col}')

    include_columns = [str(col) for col in external_cfg.get('include_columns', []) if str(col) in raw.columns]
    if not include_columns:
        excluded = {id_col, *[str(col) for col in external_cfg.get('exclude_columns', [])]}
        include_columns = [col for col in raw.columns if col not in excluded]

    if not include_columns:
        raise ValueError('No usable columns were selected from external_patient_features.path')

    prefix = str(external_cfg.get('prefix', ''))
    rename_columns = {str(k): str(v) for k, v in (external_cfg.get('rename_columns', {}) or {}).items()}

    feature_df = raw[[id_col] + include_columns].copy()
    feature_df[id_col] = feature_df[id_col].map(normalize_patient_id)
    feature_df = feature_df.rename(columns={id_col: 'patient_id'})

    final_feature_columns = []
    for column in include_columns:
        out_column = rename_columns.get(column, column)
        if prefix:
            out_column = f'{prefix}{out_column}'
        if out_column != column:
            feature_df[out_column] = feature_df[column]
            feature_df = feature_df.drop(columns=[column])
        final_feature_columns.append(out_column)

    for column in final_feature_columns:
        feature_df[column] = pd.to_numeric(feature_df[column], errors='coerce')

    agg_map = {column: first_non_null for column in final_feature_columns}
    feature_df = feature_df.groupby('patient_id', as_index=False).agg(agg_map)

    merged = base.merge(feature_df, on='patient_id', how='left')
    manifest = {
        'enabled': True,
        'path': path,
        'feature_columns': final_feature_columns,
        'n_features': len(final_feature_columns),
        'n_patients_with_any_external_feature': int(merged[final_feature_columns].notna().any(axis=1).sum()) if final_feature_columns else 0,
    }
    return merged.drop(columns=['patient_id']), manifest
