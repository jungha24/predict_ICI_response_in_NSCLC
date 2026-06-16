from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


def _fit_patient_clustering(
    patient_df: pd.DataFrame,
    feature_cols: list[str],
    cfg: dict,
    compute_umap: bool = True,
    fixed_k: int | None = None,
):
    pc_cfg = cfg["patient_clustering"]
    n_pcs_patient = int(pc_cfg["n_pcs_patient"])
    k_range = list(pc_cfg["k_range"])
    umap_n_neighbors = int(pc_cfg.get("umap_n_neighbors", 15))
    umap_min_dist = float(pc_cfg.get("umap_min_dist", 0.3))
    random_state = int(pc_cfg.get("random_state", cfg.get("seed", 123)))
    hierarchical_method = pc_cfg.get("hierarchical_method", "ward").lower()
    n_patients = int(patient_df.shape[0])

    if n_patients == 0:
        raise ValueError("Patient clustering requires at least one patient.")
    if len(feature_cols) == 0:
        raise ValueError("Patient clustering requires at least one selected feature.")

    x = patient_df[feature_cols].values
    x = StandardScaler().fit_transform(x)

    n_comp = min(n_pcs_patient, x.shape[0], x.shape[1])
    if n_comp < 1:
        raise ValueError("Patient clustering could not compute PCA because the feature matrix is empty.")

    pca = PCA(n_components=n_comp, random_state=random_state)
    pcs = pca.fit_transform(x)

    clustering_notes = []

    if compute_umap and n_patients >= 3:
        import umap

        reducer = umap.UMAP(
            n_neighbors=min(umap_n_neighbors, n_patients - 1),
            min_dist=umap_min_dist,
            metric="euclidean",
            random_state=random_state,
        )
        umap_emb = reducer.fit_transform(pcs)
    elif compute_umap:
        clustering_notes.append("UMAP skipped because fewer than 3 patients were available.")
        umap_emb = np.zeros((n_patients, 2), dtype=float)
        umap_emb[:, 0] = pcs[:, 0]
        if pcs.shape[1] > 1:
            umap_emb[:, 1] = pcs[:, 1]
    else:
        umap_emb = None

    if hierarchical_method in {"ward", "ward.d2"}:
        linkage_method = "ward"
    else:
        linkage_method = hierarchical_method

    silhouette_rows = []
    best_k = 1
    best_sil = np.nan

    if n_patients == 1:
        clustering_notes.append("Hierarchical clustering reduced to a single cluster because only one patient was available.")
        cluster_labels = np.array([1], dtype=int)
    else:
        Z = linkage(x, method=linkage_method, metric="euclidean")
        if fixed_k is None:
            valid_k_range = sorted({int(k) for k in k_range if 2 <= int(k) <= n_patients})
            silhouette_k_range = [k for k in valid_k_range if k < n_patients]

            for k in valid_k_range:
                labels = fcluster(Z, t=k, criterion="maxclust")
                if k in silhouette_k_range and len(np.unique(labels)) >= 2:
                    sil = silhouette_score(x, labels)
                else:
                    sil = np.nan

                silhouette_rows.append({"k": k, "avg_silhouette": float(sil) if pd.notna(sil) else np.nan})

                if pd.notna(sil) and (pd.isna(best_sil) or sil > best_sil):
                    best_sil = float(sil)
                    best_k = k

            if not silhouette_k_range:
                clustering_notes.append("Silhouette-based k selection skipped because fewer than 3 patients were available.")
            elif pd.isna(best_sil):
                clustering_notes.append("Silhouette scores were unavailable for the configured k_range; using a fallback cluster count.")

            if pd.isna(best_sil):
                if valid_k_range:
                    best_k = min(valid_k_range[0], n_patients)
                else:
                    best_k = min(max(k_range[0], 1), n_patients)
        else:
            best_k = int(fixed_k)
            if not 1 <= best_k <= n_patients:
                raise ValueError(f"fixed_k must be between 1 and {n_patients}, got {best_k}")
            clustering_notes.append(f"Cluster count fixed at k={best_k}.")

        cluster_labels = fcluster(Z, t=best_k, criterion="maxclust")

        if 2 <= best_k < n_patients and len(np.unique(cluster_labels)) >= 2:
            best_sil = float(silhouette_score(x, cluster_labels))
        else:
            best_sil = np.nan

        silhouette_rows.append({"k": int(best_k), "avg_silhouette": float(best_sil) if pd.notna(best_sil) else np.nan})

    metrics = {
        "n_patients": n_patients,
        "n_features": int(len(feature_cols)),
        "best_k": int(best_k),
        "best_silhouette": float(best_sil) if pd.notna(best_sil) else np.nan,
        "selected_feature_names": feature_cols,
        "clustering_notes": clustering_notes,
    }

    return pcs, pca.explained_variance_ratio_, umap_emb, cluster_labels.astype(int), pd.DataFrame(silhouette_rows), metrics


def run_patient_clustering(patient_df: pd.DataFrame, feature_cols: list[str], cfg: dict):
    pcs, explained_variance_ratio, umap_emb, cluster_labels, sil_df, metrics = _fit_patient_clustering(
        patient_df,
        feature_cols,
        cfg,
        compute_umap=True,
    )

    patient_pca_variance_df = pd.DataFrame({
        "pc": [f"PC{i+1}" for i in range(pcs.shape[1])],
        "variance_ratio": explained_variance_ratio,
    })
    patient_pca_variance_df["cumulative_variance_ratio"] = patient_pca_variance_df["variance_ratio"].cumsum()

    pca_df = pd.DataFrame(
        pcs[:, : min(10, pcs.shape[1])],
        columns=[f"PC{i+1}" for i in range(min(10, pcs.shape[1]))]
    )
    pca_df["analysis_id"] = patient_df["analysis_id"].values

    umap_df = pd.DataFrame(umap_emb, columns=["UMAP1", "UMAP2"])
    umap_df["analysis_id"] = patient_df["analysis_id"].values

    cluster_df = pd.DataFrame({
        "analysis_id": patient_df["analysis_id"].values,
        "patient_cluster": cluster_labels,
    })

    return pca_df, umap_df, cluster_df, sil_df, metrics, patient_pca_variance_df


def score_patient_clustering(patient_df: pd.DataFrame, feature_cols: list[str], cfg: dict):
    _, _, _, cluster_labels, sil_df, metrics = _fit_patient_clustering(
        patient_df,
        feature_cols,
        cfg,
        compute_umap=False,
        fixed_k=None,
    )

    cluster_df = pd.DataFrame({
        "analysis_id": patient_df["analysis_id"].values,
        "patient_cluster": cluster_labels,
    })

    return cluster_df, sil_df, metrics


def score_patient_clustering_fixed_k(
    patient_df: pd.DataFrame,
    feature_cols: list[str],
    cfg: dict,
    fixed_k: int,
):
    _, _, _, cluster_labels, sil_df, metrics = _fit_patient_clustering(
        patient_df,
        feature_cols,
        cfg,
        compute_umap=False,
        fixed_k=int(fixed_k),
    )

    cluster_df = pd.DataFrame({
        "analysis_id": patient_df["analysis_id"].values,
        "patient_cluster": cluster_labels,
    })

    return cluster_df, sil_df, metrics
