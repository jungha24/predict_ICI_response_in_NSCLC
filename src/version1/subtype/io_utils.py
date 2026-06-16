from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, List, Tuple

import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.io import mmread


def read_yaml(path: str | Path) -> dict:
    path = Path(path)
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def write_yaml(obj: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(obj, f, sort_keys=False, allow_unicode=True)


def read_json(path: str | Path) -> Any:
    path = Path(path)
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(obj: Any, path: str | Path, indent: int = 2) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=indent, ensure_ascii=False)


def read_mtx(path: str | Path) -> sp.csr_matrix:
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as f:
            mat = mmread(f)
    else:
        mat = mmread(path)
    return mat.tocsr()


def read_tsv_list(path: str | Path) -> List[str]:
    path = Path(path)
    return pd.read_csv(path, sep="\t", header=None)[0].astype(str).tolist()


def load_inputs(cfg: dict) -> Tuple[sp.csr_matrix, List[str], List[str], pd.DataFrame]:
    counts = read_mtx(cfg["paths"]["counts_mtx"])
    genes = read_tsv_list(cfg["paths"]["genes_tsv"])
    barcodes = read_tsv_list(cfg["paths"]["barcodes_tsv"])
    meta = pd.read_parquet(cfg["paths"]["cell_meta"])

    if counts.shape[0] != len(genes):
        raise ValueError(
            f"Gene count mismatch: matrix has {counts.shape[0]} rows, genes.tsv has {len(genes)} lines"
        )
    if counts.shape[1] != len(barcodes):
        raise ValueError(
            f"Barcode count mismatch: matrix has {counts.shape[1]} cols, barcodes.tsv has {len(barcodes)} lines"
        )
    if meta.shape[0] != len(barcodes):
        raise ValueError(
            f"cell_meta row mismatch: parquet has {meta.shape[0]} rows, barcodes.tsv has {len(barcodes)} lines"
        )

    meta = meta.copy()
    if "cell_id" not in meta.columns:
        raise ValueError("cell_meta.parquet must contain 'cell_id' column")
    meta["cell_id"] = meta["cell_id"].astype(str)
    meta = meta.set_index("cell_id").loc[barcodes].reset_index()

    return counts, genes, barcodes, meta


def ensure_outdir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_run_outputs(run_dir: str | Path) -> dict:
    run_dir = Path(run_dir)
    out = {"run_dir": run_dir}

    csv_files = {
        "patient_features": "patient_features.csv",
        "celltype_proportions_long": "celltype_proportions_long.csv",
        "valid_pseudobulk_groups": "valid_pseudobulk_groups.csv",
        "pseudobulk_group_metadata": "pseudobulk_group_metadata.csv",
        "pca_explained_variance": "pca_explained_variance.csv",
        "nmf_module_scores": "nmf_module_scores.csv",
        "nmf_gene_lists": "nmf_gene_lists.csv",
        "prop_related_df": "prop_related_df.csv",
        "resid_signal_df": "resid_signal_df.csv",
        "pca_nmf_redundancy": "pca_nmf_redundancy.csv",
        "patient_pca": "patient_pca.csv",
        "patient_pca_explained_variance": "patient_pca_explained_variance.csv",
        "patient_umap": "patient_umap.csv",
        "patient_clusters": "patient_clusters.csv",
        "silhouette_by_k": "silhouette_by_k.csv",
    }
    meta_files = {
        "config": "config.yaml",
        "selected_features": "selected_features.json",
        "nmf_top_genes": "nmf_top_genes.json",
        "cohort_info": "cohort_info.json",
        "cohort_summary": "cohort_summary.json",
        "metrics": "metrics.json",
    }

    for key, fname in csv_files.items():
        fpath = run_dir / fname
        out[key] = pd.read_csv(fpath) if fpath.exists() else None

    if out["patient_features"] is None:
        legacy = run_dir / "patient_features_raw.csv"
        out["patient_features"] = pd.read_csv(legacy) if legacy.exists() else None

    out["patient_features_raw"] = out["patient_features"]
    out["patient_features_batch_adjusted"] = out["patient_features"]

    for key, fname in meta_files.items():
        fpath = run_dir / fname
        if not fpath.exists():
            out[key] = None
        elif fpath.suffix == ".yaml":
            out[key] = read_yaml(fpath)
        else:
            out[key] = read_json(fpath)

    return out
