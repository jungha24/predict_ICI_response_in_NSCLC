"""Small dataclasses plus config-to-analysis-plan resolution helpers."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd


DEFAULT_COLUMN_MAP = {
    "patient_id": "sample_id",
    "sex": "Sex M: 1  F: 2",
    "age": "Age_IO",
    "histology": "Histology_mod",
    "smoking": "Smoking Never: 0 Ex: 1 Current: 2",
    "ecog": "ECOG _PS",
    "drug": "Drug",
    "pd_event": "PD_Event",
    "death_event": "Death_Event",
    "pfs_days": "PFS (Days)",
    "pfs_event": None,
    "recist": "RECIST response",
    "benefit": "Binarized response",
    "egfr_status": "EGFR_Mutation_Status_mod",
    "io_line": "IO_Line",
    "previous_chemo": "Previous_palliative_chemo",
    "previous_target": "Previous_palliative_target",
    "pd_l1_tps": "PD-L1_TPS_mod",
}

CLINICAL_EXCLUSION_ALIASES = {
    "drug": {"drug", "drug_class"},
}


@dataclass(frozen=True)
class FeatureBlock:
    name: str
    columns: list[str]
    source_columns: list[str]
    transform: str
    prefix: str


@dataclass(frozen=True)
class AnalysisSpec:
    block: str
    name: str
    clinical_feature_set: Optional[str]
    immune_blocks: list[str]
    description: str = ""


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


@dataclass
class PreparedDataset:
    analysis_df: pd.DataFrame
    clinical_feature_sets: dict[str, list[str]]
    immune_blocks: dict[str, FeatureBlock]
    feature_dictionary: pd.DataFrame
    carry_columns: list[str]
    diagnostic_tables: dict[str, pd.DataFrame]
    manifest: dict[str, Any]


def resolve_clinical_feature_sets(cfg: dict, available_columns: list[str]) -> dict[str, list[str]]:
    # Resolve named clinical covariate sets once so primary/secondary analyses can reuse them.
    clinical_cfg = cfg.get("clinical", {})
    feature_sets_cfg = clinical_cfg.get("feature_sets", {})
    if not feature_sets_cfg:
        default_covariates = clinical_cfg.get("modeling_covariates", [])
        if not default_covariates:
            raise ValueError("clinical.feature_sets or clinical.modeling_covariates must be provided.")
        feature_sets_cfg = {"default": default_covariates}

    excluded_cfg = clinical_cfg.get("exclude_from_modeling", []) or []
    if isinstance(excluded_cfg, str):
        excluded_cfg = [excluded_cfg]
    excluded_columns: set[str] = set()
    for column in excluded_cfg:
        column_name = str(column)
        excluded_columns.add(column_name)
        excluded_columns.update(CLINICAL_EXCLUSION_ALIASES.get(column_name, set()))

    resolved: dict[str, list[str]] = {}
    missing_messages = []
    for set_name, columns in feature_sets_cfg.items():
        filtered_columns = [col for col in columns if col not in excluded_columns]
        missing = [col for col in filtered_columns if col not in available_columns]
        if missing:
            missing_messages.append(f"{set_name}: {missing}")
        resolved[set_name] = [col for col in filtered_columns if col in available_columns]

    if missing_messages:
        raise ValueError(
            "Some configured clinical feature sets reference unavailable columns:\n"
            + "\n".join(missing_messages)
        )
    return resolved


def build_analysis_specs(cfg: dict) -> dict[str, list[AnalysisSpec]]:
    # Convert the YAML analysis plan into lightweight specs the pipeline can iterate over.
    analyses_cfg = cfg.get("analyses", {})
    if not analyses_cfg:
        raise ValueError("analyses must be defined in the modeling config.")

    out: dict[str, list[AnalysisSpec]] = {}
    for block_name, specs in analyses_cfg.items():
        block_specs: list[AnalysisSpec] = []
        for spec in specs:
            if "name" not in spec:
                raise ValueError(f"Each analysis entry in '{block_name}' must define a name.")
            immune_blocks = spec.get("immune_blocks", []) or []
            if isinstance(immune_blocks, str):
                immune_blocks = [immune_blocks]
            block_specs.append(
                AnalysisSpec(
                    block=block_name,
                    name=str(spec["name"]),
                    clinical_feature_set=spec.get("clinical_feature_set"),
                    immune_blocks=[str(block_name) for block_name in immune_blocks],
                    description=str(spec.get("description", "")),
                )
            )
        out[block_name] = block_specs
    return out


def _slugify(text: str) -> str:
    value = re.sub(r"[^0-9A-Za-z]+", "_", str(text).strip())
    value = re.sub(r"_+", "_", value).strip("_")
    return value.lower()


def resolve_endpoint_specs(
    cfg: dict,
    available_columns: list[str],
) -> tuple[list[BinaryEndpointSpec], list[SurvivalEndpointSpec]]:
    # Endpoint specs are configurable so binary benefit can be added without hard-coding new loops.
    endpoints_cfg = cfg.get("modeling", {}).get("endpoints", {})
    binary_cfg = endpoints_cfg.get(
        "binary",
        [
            {
                "name": "Binarized_response",
                "outcome_col": "benefit",
                "description": "Clinical benefit endpoint from the precomputed Binarized response column.",
            },
        ],
    )
    survival_cfg = endpoints_cfg.get(
        "survival",
        [
            {
                "name": "PFS",
                "time_col": "pfs_days",
                "event_col": "pfs_event",
                "description": "Progression-free survival with event defined as progression or death.",
            }
        ],
    )

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
