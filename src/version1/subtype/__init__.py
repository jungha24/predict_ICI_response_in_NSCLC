"""
Subtype discovery package.
Recommended reading order:
1. cohort.py
2. pseudobulk.py
3. features.py
4. feature_selection.py
5. batch_adjust.py
6. clustering.py
7. plotting.py
"""

from .io_utils import (
    read_yaml,
    write_yaml,
    read_json,
    write_json,
    read_mtx,
    read_tsv_list,
    load_inputs,
    ensure_outdir,
    load_run_outputs,
)

from .cohort import canonicalize_metadata, build_discovery_cohort, summarize_discovery_cohort
from .pseudobulk import GROUP_SEP, make_group_key, split_group_key, get_valid_pseudobulk_groups, aggregate_pseudobulk
from .features import (
    build_patient_metadata,
    make_proportion_features,
    make_pseudobulk_pca_features,
    make_pseudobulk_nmf_features,
    fill_missing_values,
)
from .feature_selection import (
    check_prop_related_features,
    check_feature_residual_signal,
    check_pca_nmf_redundancy,
    select_features_for_clustering,
)
from .batch_adjust import resolve_covariates, make_design_matrix, residualize_features, residualize_matrix_rows
from .clustering import run_patient_clustering, score_patient_clustering, score_patient_clustering_fixed_k
from .plotting import save_basic_plots
from .scanpy_pseudobulk import run_scanpy_pca_workflow, run_combat_seq_nmf_workflow

__all__ = [
    "read_yaml",
    "write_yaml",
    "read_json",
    "write_json",
    "read_mtx",
    "read_tsv_list",
    "load_inputs",
    "ensure_outdir",
    "load_run_outputs",
    "canonicalize_metadata",
    "build_discovery_cohort",
    "summarize_discovery_cohort",
    "GROUP_SEP",
    "make_group_key",
    "split_group_key",
    "get_valid_pseudobulk_groups",
    "aggregate_pseudobulk",
    "build_patient_metadata",
    "make_proportion_features",
    "make_pseudobulk_pca_features",
    "make_pseudobulk_nmf_features",
    "fill_missing_values",
    "check_prop_related_features",
    "check_feature_residual_signal",
    "check_pca_nmf_redundancy",
    "select_features_for_clustering",
    "resolve_covariates",
    "make_design_matrix",
    "residualize_features",
    "residualize_matrix_rows",
    "run_scanpy_pca_workflow",
    "run_combat_seq_nmf_workflow",
    "run_patient_clustering",
    "score_patient_clustering",
    "score_patient_clustering_fixed_k",
    "save_basic_plots",
]
