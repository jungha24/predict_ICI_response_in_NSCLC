from __future__ import annotations

import json
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import NMF

from .batch_adjust import residualize_matrix_rows
from .features import log_cpm
from .pseudobulk import split_group_key


def _import_scanpy():
    try:
        import scanpy as sc
        import anndata as ad
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "scanpy/anndata are required for the Scanpy pseudobulk workflow. "
            "Please run this pipeline in the dedicated Scanpy environment."
        ) from e
    return sc, ad


def _import_plotting():
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "matplotlib/seaborn are required to generate PCA/NMF review plots."
        ) from e
    return plt, sns


def _sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def _to_dense_float_matrix(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix.toarray().astype(float, copy=False)
    return np.asarray(matrix, dtype=float)


def _subset_pseudobulk_by_celltype(
    pb_counts: sp.spmatrix,
    pb_group_names: List[str],
    pb_meta: pd.DataFrame,
    celltype: str,
) -> Tuple[sp.csr_matrix, pd.DataFrame, List[int]]:
    cols = [i for i, g in enumerate(pb_group_names) if split_group_key(g)[1] == celltype]
    meta_ct = pb_meta.iloc[cols].reset_index(drop=True).copy()
    counts_ct = pb_counts[:, cols].tocsr()
    return counts_ct, meta_ct, cols


def _resolve_pc_selection(cfg: dict) -> Dict[str, int]:
    pca_cfg = cfg.get("pca_features", {})
    selection = {}

    selection_file = pca_cfg.get("pc_selection_file")
    if selection_file:
        with open(selection_file, "r", encoding="utf-8") as f:
            payload = json.load(f)
        if isinstance(payload, dict):
            selection.update({str(k): int(v) for k, v in payload.items()})

    inline_selection = pca_cfg.get("pc_selection", {})
    if isinstance(inline_selection, dict):
        selection.update({str(k): int(v) for k, v in inline_selection.items()})

    return selection


def _resolve_selected_nmf_modules(cfg: dict) -> List[str]:
    nmf_cfg = cfg.get("nmf_features", {})
    selected = [str(x) for x in nmf_cfg.get("selected_modules", [])]

    selected_file = nmf_cfg.get("selected_modules_file")
    if selected_file:
        with open(selected_file, "r", encoding="utf-8") as f:
            payload = json.load(f)
        if isinstance(payload, dict):
            selected.extend(str(x) for x in payload.get("selected_modules", []))
        elif isinstance(payload, list):
            selected.extend(str(x) for x in payload)

    return list(dict.fromkeys(selected))


def _resolve_batch_key_and_covariates(sample_meta: pd.DataFrame, covariates: List[str]) -> Tuple[Optional[str], List[str]]:
    present = [c for c in covariates if c in sample_meta.columns]
    if not present:
        return None, []

    if "batch" in present:
        batch_key = "batch"
    else:
        categorical = [c for c in present if not pd.api.types.is_numeric_dtype(sample_meta[c])]
        batch_key = categorical[0] if categorical else present[0]

    other = [c for c in present if c != batch_key]
    return batch_key, other


def _save_scree_plot(
    variance_df: pd.DataFrame,
    celltype: str,
    chosen_npcs: int,
    outpath: Path,
) -> None:
    plt, _ = _import_plotting()
    plt.figure(figsize=(5, 4))
    plt.plot(variance_df["pc"], variance_df["variance_ratio"], marker="o")
    plt.axvline(chosen_npcs, color="tab:red", linestyle="--", linewidth=1)
    plt.xlabel("Principal component")
    plt.ylabel("Explained variance ratio")
    plt.title(f"{celltype} scree plot")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def _save_nmf_heatmap(
    component_weights: np.ndarray,
    gene_names: List[str],
    module_names: List[str],
    celltype: str,
    outpath: Path,
    top_n: int,
) -> None:
    plt, sns = _import_plotting()

    heatmap_rows = []
    heatmap_index = []
    union_genes: List[str] = []
    for i, module_name in enumerate(module_names):
        order = np.argsort(-component_weights[i])[:top_n]
        genes = [gene_names[j] for j in order]
        union_genes.extend(genes)
        heatmap_rows.append((module_name, genes))

    union_genes = list(dict.fromkeys(union_genes))
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    heatmap = np.zeros((len(module_names), len(union_genes)), dtype=float)

    for row_idx, (module_name, genes) in enumerate(heatmap_rows):
        heatmap_index.append(module_name)
        for gene in genes:
            heatmap[row_idx, union_genes.index(gene)] = component_weights[row_idx, gene_to_idx[gene]]

    plt.figure(figsize=(min(18, max(8, len(union_genes) * 0.35)), max(4, len(module_names) * 0.7)))
    sns.heatmap(
        pd.DataFrame(heatmap, index=heatmap_index, columns=union_genes),
        cmap="viridis",
        cbar_kws={"label": "Module weight"},
    )
    plt.title(f"{celltype} NMF top-gene heatmap")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def run_scanpy_pca_workflow(
    pb_counts: sp.spmatrix,
    pb_group_names: List[str],
    pb_meta: pd.DataFrame,
    genes: List[str],
    cfg: dict,
    outdir: Path,
    covariates: List[str],
) -> Tuple[Optional[pd.DataFrame], pd.DataFrame]:
    sc, ad = _import_scanpy()

    pca_cfg = cfg.get("pca_features", {})
    pca_celltypes = cfg["celltypes"].get("pseudobulk_pca", [])
    min_patients = int(pca_cfg.get("min_patients_per_celltype_pca", 0))
    default_npcs = int(pca_cfg.get("npcs_per_celltype", 5))
    max_npcs = int(pca_cfg.get("max_npcs_to_compute", max(default_npcs, 20)))
    use_hvg = bool(pca_cfg.get("use_hvg", True))
    hvg_n_top_genes = int(pca_cfg.get("hvg_n_top_genes", 2000))
    hvg_flavor = str(pca_cfg.get("hvg_flavor", "seurat"))
    selection = _resolve_pc_selection(cfg)

    scree_dir = outdir / "pca_scree_plots"
    scree_dir.mkdir(parents=True, exist_ok=True)
    pca_dir = outdir / "pca_tables"
    pca_dir.mkdir(parents=True, exist_ok=True)

    feature_tables = []
    variance_tables = []

    for celltype in pca_celltypes:
        counts_ct, meta_ct, cols = _subset_pseudobulk_by_celltype(pb_counts, pb_group_names, pb_meta, celltype)
        if counts_ct.shape[1] < min_patients:
            continue

        meta_ct = meta_ct.copy()
        meta_ct.index = meta_ct["group_key"].astype(str)
        adata = ad.AnnData(
            X=counts_ct.T.copy(),
            obs=meta_ct,
            var=pd.DataFrame(index=np.array(genes, dtype=str)),
        )
        adata.var_names = np.array(genes, dtype=str)
        adata.layers["counts"] = adata.X.copy()

        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)

        hvg_mask = np.ones(adata.n_vars, dtype=bool)
        if use_hvg and adata.n_vars >= 2:
            n_top = min(hvg_n_top_genes, adata.n_vars)
            if n_top >= 2:
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_top,
                    flavor=hvg_flavor,
                    subset=False,
                    inplace=True,
                )
                hvg_mask = adata.var["highly_variable"].fillna(False).to_numpy(dtype=bool)
                if hvg_mask.sum() < 2:
                    hvg_mask = np.ones(adata.n_vars, dtype=bool)
        hvg_gene_table = pd.DataFrame(
            {
                "gene": adata.var_names.astype(str),
                "highly_variable": hvg_mask,
            }
        )

        log_cpm_matrix = _to_dense_float_matrix(adata.X).copy()
        adata.layers["log_cpm"] = log_cpm_matrix

        corrected = residualize_matrix_rows(
            log_cpm_matrix.T,
            meta_ct,
            covariates,
        ).T if covariates else log_cpm_matrix.copy()
        adata.layers["log_cpm_corrected"] = corrected.copy()
        adata.X = corrected

        gene_var = np.var(_to_dense_float_matrix(adata.X), axis=0)
        keep = gene_var > 0
        if use_hvg:
            keep = keep & hvg_mask
        adata = adata[:, keep].copy()
        if adata.n_vars < 2:
            continue

        n_comps = min(max_npcs, adata.n_obs, adata.n_vars)
        if n_comps < 1:
            continue

        sc.pp.scale(adata)
        sc.pp.pca(adata, n_comps=n_comps)

        variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=float)
        variance_df = pd.DataFrame({
            "celltype": celltype,
            "pc": np.arange(1, variance_ratio.size + 1),
            "variance_ratio": variance_ratio,
            "cumulative_variance_ratio": np.cumsum(variance_ratio),
        })
        variance_tables.append(variance_df)

        chosen_npcs = min(int(selection.get(celltype, default_npcs)), variance_ratio.size)
        chosen_npcs = max(chosen_npcs, 1)
        scores = np.asarray(adata.obsm["X_pca"][:, :chosen_npcs], dtype=float)

        feature_df = pd.DataFrame(
            scores,
            columns=[f"{celltype}_PC{i+1}" for i in range(scores.shape[1])],
        )
        feature_df["analysis_id"] = meta_ct["analysis_id"].values
        feature_tables.append(feature_df)

        safe_ct = _sanitize_name(celltype)
        hvg_gene_table.to_csv(pca_dir / f"{safe_ct}_hvg_genes.csv", index=False)
        variance_df.to_csv(pca_dir / f"{safe_ct}_pca_variance.csv", index=False)
        pd.DataFrame(
            np.asarray(adata.obsm["X_pca"], dtype=float),
            columns=[f"PC{i+1}" for i in range(np.asarray(adata.obsm["X_pca"]).shape[1])],
        ).assign(
            analysis_id=meta_ct["analysis_id"].values,
            group_key=meta_ct["group_key"].values,
        ).to_csv(pca_dir / f"{safe_ct}_pca_scores.csv", index=False)
        _save_scree_plot(
            variance_df,
            celltype=celltype,
            chosen_npcs=chosen_npcs,
            outpath=scree_dir / f"{safe_ct}_scree.png",
        )

    merged = None
    if feature_tables:
        merged = feature_tables[0]
        for df in feature_tables[1:]:
            merged = merged.merge(df, on="analysis_id", how="outer")

    variance_out = pd.concat(variance_tables, ignore_index=True) if variance_tables else pd.DataFrame(
        columns=["celltype", "pc", "variance_ratio", "cumulative_variance_ratio"]
    )
    return merged, variance_out


def _assess_combat_seq_feasibility(meta_df: pd.DataFrame, batch_key: Optional[str]) -> Tuple[bool, Optional[str]]:
    if not batch_key or batch_key not in meta_df.columns:
        return False, "no batch key configured"

    batch = meta_df[batch_key].fillna("NA").astype(str)
    n_unique = int(batch.nunique())
    if n_unique < 2:
        return False, "fewer than two batches are present"

    counts = batch.value_counts(dropna=False)
    singleton_batches = counts[counts < 2]
    if not singleton_batches.empty:
        singleton_label = ", ".join(f"{idx}:{int(val)}" for idx, val in singleton_batches.items())
        return False, f"singleton batches present ({singleton_label})"

    return True, None


def _run_combat_seq(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    batch_key: str,
    covariates: List[str],
    script_path: Path,
) -> pd.DataFrame:
    rscript = "Rscript"
    with tempfile.TemporaryDirectory(prefix="combat_seq_") as tmpdir:
        tmpdir = Path(tmpdir)
        counts_path = tmpdir / "counts.csv"
        meta_path = tmpdir / "meta.csv"
        out_path = tmpdir / "adjusted_counts.csv"

        counts_df.to_csv(counts_path)
        meta_df.to_csv(meta_path, index=False)

        cmd = [
            rscript,
            str(script_path),
            str(counts_path),
            str(meta_path),
            str(out_path),
            batch_key,
            ",".join(covariates),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                "ComBat-seq failed.\n"
                f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
            )

        adjusted = pd.read_csv(out_path, index_col=0)

    adjusted = adjusted.loc[counts_df.index, counts_df.columns]
    return adjusted


def run_combat_seq_nmf_workflow(
    pb_counts: sp.spmatrix,
    pb_group_names: List[str],
    pb_meta: pd.DataFrame,
    genes: List[str],
    cfg: dict,
    outdir: Path,
    covariates: List[str],
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Dict[str, dict], pd.DataFrame]:
    nmf_cfg = cfg.get("nmf_features", {})
    if not nmf_cfg.get("enabled", True):
        return None, None, {}, pd.DataFrame(columns=["celltype", "module", "rank", "gene", "weight"])

    nmf_celltypes = cfg["celltypes"].get("pseudobulk_nmf", [])
    min_patients = int(nmf_cfg.get("min_patients_per_celltype_nmf", 0))
    min_genes = int(nmf_cfg.get("min_genes_nmf", 200))
    top_var_genes = int(nmf_cfg.get("nmf_top_var_genes", 1500))
    rank = int(nmf_cfg.get("rank", 3))
    top_genes_per_module = int(nmf_cfg.get("top_genes_per_module", 30))
    heatmap_top_genes = int(nmf_cfg.get("heatmap_top_genes", 15))
    selected_modules = set(_resolve_selected_nmf_modules(cfg))

    nmf_dir = outdir / "nmf_review"
    nmf_dir.mkdir(parents=True, exist_ok=True)
    heatmap_dir = nmf_dir / "heatmaps"
    heatmap_dir.mkdir(parents=True, exist_ok=True)
    counts_dir = nmf_dir / "combat_seq_counts"
    counts_dir.mkdir(parents=True, exist_ok=True)

    script_path = Path(__file__).resolve().parents[2] / "Rscripts" / "combat_seq_pseudobulk.R"
    batch_key, extra_covariates = _resolve_batch_key_and_covariates(pb_meta, covariates)

    all_module_tables = []
    selected_module_tables = []
    gene_rows = []
    nmf_meta: Dict[str, dict] = {}

    for celltype in nmf_celltypes:
        counts_ct, meta_ct, cols = _subset_pseudobulk_by_celltype(pb_counts, pb_group_names, pb_meta, celltype)
        if counts_ct.shape[1] < min_patients:
            continue

        count_df = pd.DataFrame(
            np.asarray(counts_ct.toarray(), dtype=float),
            index=np.array(genes, dtype=str),
            columns=meta_ct["group_key"].astype(str).tolist(),
        )

        combat_seq_applied = False
        combat_seq_skip_reason = None
        can_run_combat, infeasible_reason = _assess_combat_seq_feasibility(meta_ct, batch_key)

        if batch_key and can_run_combat:
            try:
                adjusted_counts = _run_combat_seq(
                    count_df,
                    meta_ct,
                    batch_key=batch_key,
                    covariates=extra_covariates,
                    script_path=script_path,
                )
                combat_seq_applied = True
            except RuntimeError as exc:
                error_text = str(exc)
                if "1 sample per batch" in error_text:
                    adjusted_counts = count_df.copy()
                    combat_seq_skip_reason = "ComBat-seq rejected singleton batches at runtime"
                else:
                    raise
        else:
            adjusted_counts = count_df.copy()
            if batch_key:
                combat_seq_skip_reason = infeasible_reason
            else:
                combat_seq_skip_reason = "batch correction disabled"

        if combat_seq_skip_reason:
            print(f"[subtype.nmf] Skipping ComBat-seq for {celltype}: {combat_seq_skip_reason}", flush=True)

        safe_ct = _sanitize_name(celltype)
        adjusted_counts.to_csv(counts_dir / f"{safe_ct}_combat_seq_counts.csv")

        adjusted_sparse = sp.csr_matrix(adjusted_counts.values)
        expr_log = log_cpm(adjusted_sparse)
        gene_var = expr_log.var(axis=1)
        keep_var = np.argsort(-gene_var)[: min(top_var_genes, len(gene_var))]
        expr_nmf = np.clip(adjusted_counts.values[keep_var, :], a_min=0.0, a_max=None).T

        if expr_nmf.shape[1] < min_genes:
            continue

        n_comp = min(rank, expr_nmf.shape[0], expr_nmf.shape[1])
        if n_comp < 1:
            continue

        model = NMF(
            n_components=n_comp,
            init="nndsvda",
            random_state=cfg.get("seed", 123),
            max_iter=1000,
        )
        W = model.fit_transform(expr_nmf)
        H = model.components_

        module_names = [f"{celltype}_NMF{i+1}" for i in range(W.shape[1])]
        module_df = pd.DataFrame(W, columns=module_names)
        module_df["analysis_id"] = meta_ct["analysis_id"].values
        module_df["group_key"] = meta_ct["group_key"].values
        all_module_tables.append(module_df)

        kept_gene_names = adjusted_counts.index.to_numpy()[keep_var].tolist()
        top_genes = {}
        for i, module_name in enumerate(module_names):
            order = np.argsort(-H[i])[:top_genes_per_module]
            genes_for_module = [kept_gene_names[j] for j in order]
            weights_for_module = [float(H[i, j]) for j in order]
            top_genes[module_name] = genes_for_module
            for rank_idx, (gene, weight) in enumerate(zip(genes_for_module, weights_for_module), start=1):
                gene_rows.append({
                    "celltype": celltype,
                    "module": module_name,
                    "rank": rank_idx,
                    "gene": gene,
                    "weight": weight,
                })

        _save_nmf_heatmap(
            component_weights=H,
            gene_names=kept_gene_names,
            module_names=module_names,
            celltype=celltype,
            outpath=heatmap_dir / f"{safe_ct}_nmf_heatmap.png",
            top_n=heatmap_top_genes,
        )

        nmf_meta[celltype] = {
            "rank_used": int(n_comp),
            "batch_key": batch_key,
            "covariates_used": extra_covariates,
            "combat_seq_applied": combat_seq_applied,
            "combat_seq_skip_reason": combat_seq_skip_reason,
            "top_genes": top_genes,
            "combat_seq_counts_file": str(counts_dir / f"{safe_ct}_combat_seq_counts.csv"),
        }

        keep_cols = [c for c in module_names if c in selected_modules]
        if keep_cols:
            selected_module_tables.append(module_df[["analysis_id"] + keep_cols].copy())

    all_modules_merged = None
    if all_module_tables:
        all_modules_merged = all_module_tables[0]
        for df in all_module_tables[1:]:
            all_modules_merged = all_modules_merged.merge(df, on=["analysis_id", "group_key"], how="outer")

    selected_modules_merged = None
    if selected_module_tables:
        selected_modules_merged = selected_module_tables[0]
        for df in selected_module_tables[1:]:
            selected_modules_merged = selected_modules_merged.merge(df, on="analysis_id", how="outer")

    gene_list_df = pd.DataFrame(gene_rows, columns=["celltype", "module", "rank", "gene", "weight"])
    return selected_modules_merged, all_modules_merged, nmf_meta, gene_list_df
