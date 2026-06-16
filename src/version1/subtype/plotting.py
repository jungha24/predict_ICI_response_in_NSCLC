from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def save_basic_plots(
    outdir: Path,
    pca_df: pd.DataFrame,
    umap_df: pd.DataFrame,
    patient_final: pd.DataFrame,
    patient_pca_variance_df: pd.DataFrame | None,
    cfg: dict,
) -> None:
    output_cfg = cfg.get("output", {})
    if not output_cfg.get("save_plots", True):
        return

    plots_dirname = output_cfg.get("plots_dirname", "plots")
    plots_dir = outdir / plots_dirname
    plots_dir.mkdir(parents=True, exist_ok=True)

    pca_plot_df = pca_df.merge(
        patient_final[["analysis_id", "batch", "patient_cluster"]],
        on="analysis_id",
        how="left",
    )
    umap_plot_df = umap_df.merge(
        patient_final[["analysis_id", "batch", "patient_cluster"]],
        on="analysis_id",
        how="left",
    )

    plt.figure(figsize=(5, 4))
    for s in sorted(pca_plot_df["batch"].dropna().unique()):
        sub = pca_plot_df[pca_plot_df["batch"] == s]
        plt.scatter(sub["PC1"], sub["PC2"], s=40, label=s)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Patient PCA by batch")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "pca_by_batch.png", dpi=200)
    plt.close()

    plt.figure(figsize=(5, 4))
    for s in sorted(pca_plot_df["patient_cluster"].dropna().unique()):
        sub = pca_plot_df[pca_plot_df["patient_cluster"] == s]
        plt.scatter(sub["PC1"], sub["PC2"], s=40, label=f"C{s}")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Patient PCA by cluster")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "pca_by_cluster.png", dpi=200)
    plt.close()

    plt.figure(figsize=(5, 4))
    for s in sorted(umap_plot_df["batch"].dropna().unique()):
        sub = umap_plot_df[umap_plot_df["batch"] == s]
        plt.scatter(sub["UMAP1"], sub["UMAP2"], s=40, label=s)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Patient UMAP by batch")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "umap_by_batch.png", dpi=200)
    plt.close()

    plt.figure(figsize=(5, 4))
    for s in sorted(umap_plot_df["patient_cluster"].dropna().unique()):
        sub = umap_plot_df[umap_plot_df["patient_cluster"] == s]
        plt.scatter(sub["UMAP1"], sub["UMAP2"], s=40, label=f"C{s}")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Patient UMAP by cluster")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "umap_by_cluster.png", dpi=200)
    plt.close()

    if patient_pca_variance_df is not None and not patient_pca_variance_df.empty:
        chosen_npcs = int(cfg.get("patient_clustering", {}).get("n_pcs_patient", patient_pca_variance_df.shape[0]))
        x = range(1, patient_pca_variance_df.shape[0] + 1)

        plt.figure(figsize=(5, 4))
        plt.plot(x, patient_pca_variance_df["variance_ratio"], marker="o")
        plt.axvline(min(chosen_npcs, patient_pca_variance_df.shape[0]), color="tab:red", linestyle="--", linewidth=1)
        plt.xlabel("Principal component")
        plt.ylabel("Explained variance ratio")
        plt.title("Patient PCA scree plot")
        plt.tight_layout()
        plt.savefig(plots_dir / "patient_pca_scree.png", dpi=200)
        plt.close()
