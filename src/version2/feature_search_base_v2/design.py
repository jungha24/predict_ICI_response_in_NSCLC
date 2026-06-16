from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd


@dataclass(frozen=True)
class BinaryEndpointSpec:
    name: str
    outcome_col: str
    slug: str
    description: str = ""


@dataclass(frozen=True)
class SurvivalEndpointSpec:
    name: str
    time_col: str
    event_col: str
    slug: str
    description: str = ""


@dataclass(frozen=True)
class SearchDataset:
    analysis_df: pd.DataFrame
    clinical_feature_sets: dict[str, list[str]]
    candidate_features: list[str]
    carry_columns: list[str]
    feature_catalog: pd.DataFrame
    manifest: dict[str, Any]


def _slugify(text: str) -> str:
    import re
    value = re.sub(r"[^0-9A-Za-z]+", "_", str(text).strip())
    value = re.sub(r"_+", "_", value).strip("_")
    return value.lower()


def resolve_clinical_feature_sets(cfg: dict, available_columns: list[str]) -> dict[str, list[str]]:
    clinical_cfg = cfg.get("clinical", {})
    feature_sets_cfg = clinical_cfg.get("feature_sets", {})
    if not feature_sets_cfg:
        raise ValueError("clinical.feature_sets must be provided.")

    excluded_cfg = clinical_cfg.get("exclude_from_modeling", []) or []
    if isinstance(excluded_cfg, str):
        excluded_cfg = [excluded_cfg]
    excluded = {str(x) for x in excluded_cfg}
    if "drug" in excluded:
        excluded.add("drug_class")

    resolved: dict[str, list[str]] = {}
    missing_messages = []
    for set_name, columns in feature_sets_cfg.items():
        filtered = [col for col in columns if col not in excluded]
        missing = [col for col in filtered if col not in available_columns]
        if missing:
            missing_messages.append(f"{set_name}: {missing}")
        resolved[set_name] = [col for col in filtered if col in available_columns]

    if missing_messages:
        raise ValueError(
            "Some configured clinical feature sets reference unavailable columns:\n" + "\n".join(missing_messages)
        )
    return resolved


def resolve_endpoint_specs(cfg: dict, available_columns: list[str]) -> tuple[list[BinaryEndpointSpec], list[SurvivalEndpointSpec]]:
    endpoints_cfg = cfg.get("modeling", {}).get("endpoints", {})
    binary_cfg = endpoints_cfg.get("binary", [])
    survival_cfg = endpoints_cfg.get("survival", [])

    binary_specs: list[BinaryEndpointSpec] = []
    for spec in binary_cfg:
        outcome_col = str(spec["outcome_col"])
        if outcome_col not in available_columns:
            continue
        binary_specs.append(
            BinaryEndpointSpec(
                name=str(spec.get("name", outcome_col)),
                outcome_col=outcome_col,
                slug=str(spec.get("slug", _slugify(spec.get("name", outcome_col)))),
                description=str(spec.get("description", "")),
            )
        )

    survival_specs: list[SurvivalEndpointSpec] = []
    for spec in survival_cfg:
        time_col = str(spec["time_col"])
        event_col = str(spec["event_col"])
        if time_col not in available_columns or event_col not in available_columns:
            continue
        survival_specs.append(
            SurvivalEndpointSpec(
                name=str(spec.get("name", time_col)),
                time_col=time_col,
                event_col=event_col,
                slug=str(spec.get("slug", _slugify(spec.get("name", time_col)))),
                description=str(spec.get("description", "")),
            )
        )

    return binary_specs, survival_specs
