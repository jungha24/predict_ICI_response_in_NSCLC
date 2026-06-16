### 20260311
### baseline subtype discovery
### steps required fine-tuning "#<<<"

### leveraged nsclc n73 in baseline and 10 heatlhy samples

### [1] leave patients/celltypes of interest

keep_celltypes <- c("CD4T", "CD8T", "NK", "B", "Mono", "DC")

seu_use <- subset(
  seu_raw_merged,
  subset = !is.na(celltype) & celltype %in% keep_celltypes
)

### [2] pseudobulk: pateint id x celltype
pb_list <- AggregateExpression(
  seu_use,
  group.by = c("patient_id", "celltype"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

pb_counts <- pb_list$RNA
dim(pb_counts)
head(colnames(pb_counts))

### [3] calculate cell type proportion 
cell_meta <- seu_use@meta.data %>%
  dplyr::select(patient_id, celltype)

prop_df <- cell_meta %>%
  dplyr::count(patient_id, celltype) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

prop_wide <- prop_df %>%
  tidyr::pivot_wider(
    id_cols = patient_id,
    names_from = celltype,
    values_from = prop,
    values_fill = 0
  )

colnames(prop_wide)[-1] <- paste0("prop_", colnames(prop_wide)[-1])

### [4] make pseudobulk PCA features per cell type
make_pseudobulk_pcs <- function(pb_counts, celltype_name, npcs = 5, min_patients = 5) {
  cols <- grep(paste0("_", celltype_name, "$"), colnames(pb_counts), value = TRUE)
  
  if (length(cols) < min_patients) {
    message("Skipping ", celltype_name, ": too few patient-celltype samples")
    return(NULL)
  }
  
  mat <- pb_counts[, cols, drop = FALSE]
  patient_ids <- sub(paste0("_", celltype_name, "$"), "", colnames(mat))
  colnames(mat) <- patient_ids
  
  # edgeR normalization
  dge <- DGEList(counts = mat)
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(dge) < 50) {
    message("Skipping ", celltype_name, ": too few genes after filtering")
    return(NULL)
  }
  
  dge <- calcNormFactors(dge)
  expr <- cpm(dge, log = TRUE, prior.count = 1)
  
  # PCA on patients
  pca <- prcomp(t(expr), center = TRUE, scale. = TRUE)
  pcs <- as.data.frame(pca$x[, 1:min(npcs, ncol(pca$x)), drop = FALSE])
  pcs$patient_id <- rownames(pcs)
  
  colnames(pcs)[1:(ncol(pcs)-1)] <- paste0(celltype_name, "_PC", seq_len(ncol(pcs)-1))
  return(pcs)
}
# run per celltype
celltypes <- keep_celltypes
pc_list <- lapply(celltypes, function(ct) {
  make_pseudobulk_pcs(pb_counts = pb_counts, celltype_name = ct, npcs = 5)
})
names(pc_list) <- celltypes
pc_list <- pc_list[!sapply(pc_list, is.null)]
# merge
feature_df <- Reduce(function(x, y) full_join(x, y, by = "patient_id"), pc_list)
feature_df <- full_join(feature_df, prop_wide, by = "patient_id")
# deal with na values
impute_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}
feature_mat <- feature_df
for (j in seq_len(ncol(feature_mat))) {
  if (colnames(feature_mat)[j] != "patient_id") {
    feature_mat[[j]] <- impute_median(feature_mat[[j]])
  }
}

### [5] nmf module as feature
# extract NMF module activity from pseudobulk matrix
library(edgeR)
library(NMF)

make_pseudobulk_nmf <- function(pb_counts, celltype_name, rank = 3, min_patients = 8, min_genes = 200) {
  cols <- grep(paste0("_", celltype_name, "$"), colnames(pb_counts), value = TRUE)
  
  if (length(cols) < min_patients) {
    message("Skipping ", celltype_name, ": too few patient-celltype samples")
    return(NULL)
  }
  
  mat <- pb_counts[, cols, drop = FALSE]
  patient_ids <- sub(paste0("_", celltype_name, "$"), "", colnames(mat))
  colnames(mat) <- patient_ids
  
  # edgeR normalization
  dge <- DGEList(counts = mat)
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(dge) < min_genes) {
    message("Skipping ", celltype_name, ": too few genes after filtering")
    return(NULL)
  }
  
  dge <- calcNormFactors(dge)
  
  # NMF는 non-negative input 필요
  expr <- cpm(dge, log = FALSE, prior.count = 1)
  expr <- expr + 1e-6
  
  # NMF
  fit <- nmf(expr, rank = rank, nrun = 30, seed = 123)
  
  W <- basis(fit)   # genes x modules
  H <- coef(fit)    # modules x patients
  
  nmf_df <- as.data.frame(t(H))
  nmf_df$patient_id <- rownames(nmf_df)
  
  module_cols <- setdiff(colnames(nmf_df), "patient_id")
  colnames(nmf_df)[match(module_cols, colnames(nmf_df))] <- 
    paste0(celltype_name, "_NMF", seq_along(module_cols))
  
  # top genes도 같이 반환
  top_genes <- lapply(seq_len(ncol(W)), function(i) {
    names(sort(W[, i], decreasing = TRUE))[1:30]
  })
  names(top_genes) <- paste0(celltype_name, "_NMF", seq_len(ncol(W)))
  
  return(list(activity = nmf_df, basis = W, coef = H, top_genes = top_genes, fit = fit))
}
# run NMF per celltype
celltypes <- c("CD4T", "CD8T", "NK", "B", "Mono", "DC")

nmf_list <- lapply(celltypes, function(ct) {
  make_pseudobulk_nmf(
    pb_counts = pb_counts,
    celltype_name = ct,
    rank = 3
  )
})
names(nmf_list) <- celltypes
#extract activity and merge to patient feature
nmf_activity_list <- lapply(nmf_list, function(x) {
  if (is.null(x)) return(NULL)
  x$activity
})
nmf_activity_list <- nmf_activity_list[!sapply(nmf_activity_list, is.null)]
nmf_feature_df <- Reduce(function(x, y) full_join(x, y, by = "patient_id"), nmf_activity_list)

feature_df_all <- full_join(feature_df, nmf_feature_df, by = "patient_id")

### [6]: feature selection
### [6-1]: correlation across features
### [6-2]: comparison of proportion vs PCA, proportion vs NMF, PCA vs NMF within same cell type
### [6-3]: find PC/NMF explained by proportion
### [6-4]: check overlap between PCA and NMF
### [6-5]: remove/residualize overlaping features 
### [6-6]: finalized feature set
# feature type 분류 함수
classify_feature_types <- function(feature_names) {
  data.frame(
    feature = feature_names,
    feature_type = dplyr::case_when(
      grepl("^prop_", feature_names) ~ "proportion",
      grepl("_PC[0-9]+$", feature_names) ~ "pca",
      grepl("_NMF[0-9]+$", feature_names) ~ "nmf",
      grepl("_NMF[0-9]+_residProp$", feature_names) ~ "nmf_resid",
      TRUE ~ "other"
    ),
    stringsAsFactors = FALSE
  )
}
# cell type 이름 추출 함수
extract_celltype_name <- function(feature_name) {
  if (grepl("^prop_", feature_name)) {
    return(sub("^prop_", "", feature_name))
  }
  if (grepl("_PC[0-9]+$", feature_name)) {
    return(sub("_PC[0-9]+$", "", feature_name))
  }
  if (grepl("_NMF[0-9]+_residProp$", feature_name)) {
    return(sub("_NMF[0-9]+_residProp$", "", feature_name))
  }
  if (grepl("_NMF[0-9]+$", feature_name)) {
    return(sub("_NMF[0-9]+$", "", feature_name))
  }
  return(NA_character_)
}
# 6-1
get_numeric_feature_matrix <- function(df, exclude_cols = c("patient_id", "source", "collection_batch")) {
  feature_cols <- setdiff(colnames(df), exclude_cols)
  num_df <- df[, feature_cols, drop = FALSE]
  num_df <- num_df[, sapply(num_df, is.numeric), drop = FALSE]
  return(num_df)
}

num_mat <- get_numeric_feature_matrix(patient_features)
cor_mat <- cor(num_mat, use = "pairwise.complete.obs")

pheatmap::pheatmap(cor_mat, main = "Feature correlation matrix")
# summerize feature pair redundancy table
get_high_correlation_pairs <- function(cor_mat, cutoff = 0.8) {
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_df) <- c("feature1", "feature2", "correlation")
  
  cor_df <- cor_df[cor_df$feature1 != cor_df$feature2, ]
  cor_df$abs_correlation <- abs(cor_df$correlation)
  
  # duplicate 제거
  cor_df$key <- apply(cor_df[, c("feature1", "feature2")], 1, function(x) {
    paste(sort(x), collapse = "__")
  })
  cor_df <- cor_df[!duplicated(cor_df$key), ]
  cor_df$key <- NULL
  
  cor_df <- cor_df[order(-cor_df$abs_correlation), ]
  cor_df <- subset(cor_df, abs_correlation >= cutoff)
  rownames(cor_df) <- NULL
  cor_df
}

high_corr_all <- get_high_correlation_pairs(cor_mat, cutoff = 0.8)
high_corr_all
feature_info <- classify_feature_types(colnames(num_mat))
feature_info$celltype <- sapply(feature_info$feature, extract_celltype_name)

annotate_feature_pairs <- function(pair_df, feature_info) {
  fi1 <- feature_info
  fi2 <- feature_info
  colnames(fi1) <- c("feature1", "type1", "celltype1")
  colnames(fi2) <- c("feature2", "type2", "celltype2")
  
  out <- dplyr::left_join(pair_df, fi1, by = "feature1")
  out <- dplyr::left_join(out, fi2, by = "feature2")
  out
}

high_corr_annot <- annotate_feature_pairs(high_corr_all, feature_info)
high_corr_annot

# 6-3 # correlation
check_prop_related_features <- function(df) {
  feature_names <- colnames(df)
  prop_cols <- grep("^prop_", feature_names, value = TRUE)
  pca_cols  <- grep("_PC[0-9]+$", feature_names, value = TRUE)
  nmf_cols  <- grep("_NMF[0-9]+$", feature_names, value = TRUE)
  
  out_list <- list()
  
  for (pcol in prop_cols) {
    celltype <- sub("^prop_", "", pcol)
    
    matched_pca <- grep(paste0("^", celltype, "_PC"), pca_cols, value = TRUE)
    matched_nmf <- grep(paste0("^", celltype, "_NMF"), nmf_cols, value = TRUE)
    matched_all <- c(matched_pca, matched_nmf)
    
    if (length(matched_all) == 0) next
    
    for (fcol in matched_all) {
      r <- suppressWarnings(cor(df[[pcol]], df[[fcol]], use = "pairwise.complete.obs"))
      out_list[[paste(pcol, fcol, sep = "__")]] <- data.frame(
        prop_feature = pcol,
        feature = fcol,
        celltype = celltype,
        feature_type = ifelse(grepl("_PC", fcol), "pca", "nmf"),
        correlation = r,
        abs_correlation = abs(r)
      )
    }
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out <- out[order(-out$abs_correlation), ]
  out
}

prop_related_df <- check_prop_related_features(patient_features)
prop_related_df
subset(prop_related_df, abs_correlation >= 0.7)

# 6-3 #regression
check_feature_residual_signal <- function(df) {
  feature_names <- colnames(df)
  prop_cols <- grep("^prop_", feature_names, value = TRUE)
  pca_cols  <- grep("_PC[0-9]+$", feature_names, value = TRUE)
  nmf_cols  <- grep("_NMF[0-9]+$", feature_names, value = TRUE)
  
  out_list <- list()
  
  for (pcol in prop_cols) {
    celltype <- sub("^prop_", "", pcol)
    
    matched_pca <- grep(paste0("^", celltype, "_PC"), pca_cols, value = TRUE)
    matched_nmf <- grep(paste0("^", celltype, "_NMF"), nmf_cols, value = TRUE)
    matched_all <- c(matched_pca, matched_nmf)
    
    if (length(matched_all) == 0) next
    
    for (fcol in matched_all) {
      fit <- lm(as.formula(paste0("`", fcol, "` ~ `", pcol, "`")), data = df)
      sm <- summary(fit)
      
      residual_sd <- sd(residuals(fit), na.rm = TRUE)
      raw_sd <- sd(df[[fcol]], na.rm = TRUE)
      
      out_list[[paste(pcol, fcol, sep = "__")]] <- data.frame(
        prop_feature = pcol,
        feature = fcol,
        celltype = celltype,
        feature_type = ifelse(grepl("_PC", fcol), "pca", "nmf"),
        r_squared = sm$r.squared,
        adj_r_squared = sm$adj.r.squared,
        prop_pvalue = coef(sm)[2, 4],
        raw_sd = raw_sd,
        residual_sd = residual_sd,
        residual_ratio = residual_sd / raw_sd
      )
    }
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out <- out[order(-out$r_squared), ]
  out
}

resid_signal_df <- check_feature_residual_signal(patient_features)
resid_signal_df

subset(resid_signal_df, r_squared >= 0.5 | residual_ratio <= 0.7)

# 6-4
check_pca_nmf_redundancy <- function(df) {
  feature_names <- colnames(df)
  pca_cols  <- grep("_PC[0-9]+$", feature_names, value = TRUE)
  nmf_cols  <- grep("_NMF[0-9]+$", feature_names, value = TRUE)
  
  out_list <- list()
  
  pca_info <- data.frame(feature = pca_cols, celltype = sapply(pca_cols, extract_celltype_name))
  nmf_info <- data.frame(feature = nmf_cols, celltype = sapply(nmf_cols, extract_celltype_name))
  
  for (ct in intersect(unique(pca_info$celltype), unique(nmf_info$celltype))) {
    pca_ct <- subset(pca_info, celltype == ct)$feature
    nmf_ct <- subset(nmf_info, celltype == ct)$feature
    
    for (pcol in pca_ct) {
      for (ncol in nmf_ct) {
        r <- suppressWarnings(cor(df[[pcol]], df[[ncol]], use = "pairwise.complete.obs"))
        out_list[[paste(pcol, ncol, sep = "__")]] <- data.frame(
          celltype = ct,
          pca_feature = pcol,
          nmf_feature = ncol,
          correlation = r,
          abs_correlation = abs(r)
        )
      }
    }
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out <- out[order(-out$abs_correlation), ]
  out
}

pca_nmf_red_df <- check_pca_nmf_redundancy(patient_features)
pca_nmf_red_df

subset(pca_nmf_red_df, abs_correlation >= 0.8)

# 6-5
# make feature seletion rule
# basic idea
# 1.proportion은 기본 유지
# 2.proportion에 너무 종속된 PCA/NMF는 제거 또는 residualize
# 3.PCA/NMF끼리 거의 같은 정보면 하나만 유지
# 4.NMF는 해석 가능성이 높아서, 동률이면 NMF를 남기고 PC를 빼는 전략도 가능
# 5.또는 주분석은 PCA, 해석용은 NMF를 따로 유지

## drop
select_features_for_clustering <- function(
    df,
    prop_related_df,
    resid_signal_df,
    pca_nmf_red_df,
    high_corr_cutoff = 0.8,
    prop_r2_cutoff = 0.5,
    residual_ratio_cutoff = 0.7,
    prefer = c("pca", "nmf")
) {
  prefer <- match.arg(prefer)
  
  all_features <- colnames(get_numeric_feature_matrix(df))
  
  # always keep proportions initially
  prop_feats <- grep("^prop_", all_features, value = TRUE)
  pca_feats  <- grep("_PC[0-9]+$", all_features, value = TRUE)
  nmf_feats  <- grep("_NMF[0-9]+$", all_features, value = TRUE)
  
  drop_feats <- character(0)
  
  # 1) drop PCA/NMF that are too strongly explained by same-celltype proportion
  flagged_prop <- subset(
    resid_signal_df,
    r_squared >= prop_r2_cutoff | residual_ratio <= residual_ratio_cutoff
  )
  drop_feats <- unique(c(drop_feats, flagged_prop$feature))
  
  # 2) PCA-NMF redundancy: if same biology repeated, keep only one class
  redundant_pairs <- subset(pca_nmf_red_df, abs_correlation >= high_corr_cutoff)
  
  if (nrow(redundant_pairs) > 0) {
    for (i in seq_len(nrow(redundant_pairs))) {
      pfeat <- redundant_pairs$pca_feature[i]
      nfeat <- redundant_pairs$nmf_feature[i]
      
      if (prefer == "pca") {
        drop_feats <- unique(c(drop_feats, nfeat))
      } else {
        drop_feats <- unique(c(drop_feats, pfeat))
      }
    }
  }
  
  selected <- setdiff(all_features, drop_feats)
  
  list(
    selected_features = selected,
    dropped_features = unique(drop_feats),
    kept_proportions = intersect(selected, prop_feats),
    kept_pca = intersect(selected, pca_feats),
    kept_nmf = intersect(selected, nmf_feats)
  )
}
# run
feature_selection_result <- select_features_for_clustering(
  df = patient_features,
  prop_related_df = prop_related_df,
  resid_signal_df = resid_signal_df,
  pca_nmf_red_df = pca_nmf_red_df,
  high_corr_cutoff = 0.8,
  prop_r2_cutoff = 0.5,
  residual_ratio_cutoff = 0.7,
  prefer = "pca"   # 주분석이면 PCA 우선 추천
)

feature_cols_selected <- feature_selection_result$selected_features
feature_cols_dropped <- feature_selection_result$dropped_features

## use residualized feature
make_prop_residualized_features <- function(df, resid_signal_df,
                                            r2_cutoff = 0.5,
                                            residual_ratio_cutoff = 0.7) {
  df_out <- df
  
  flagged <- subset(
    resid_signal_df,
    r_squared >= r2_cutoff | residual_ratio <= residual_ratio_cutoff
  )
  
  for (i in seq_len(nrow(flagged))) {
    pcol <- flagged$prop_feature[i]
    fcol <- flagged$feature[i]
    new_col <- paste0(fcol, "_residProp")
    
    fit <- lm(as.formula(paste0("`", fcol, "` ~ `", pcol, "`")), data = df)
    df_out[[new_col]] <- residuals(fit)
  }
  
  df_out
}

feature_df_selected <- feature_df_all[,feature_cols_selected]

### [8] patient-level metadata: merge metadata with feature matrix
patient_meta <- seu_use@meta.data %>%
  dplyr::select(patient_id, source, collection_batch) %>%
  distinct()
# check
patient_meta %>%
  dplyr::count(patient_id) %>%
  dplyr::filter(n > 1)
patient_features <- left_join(feature_df_all, patient_meta, by = "patient_id")


### [9] regress collection batch for features
# use resdiual + overall mean as new features
residualize_by_batch <- function(df, feature_cols, batch_col = "collection_batch") {
  out <- df
  
  for (feat in feature_cols) {
    form <- as.formula(paste0("`", feat, "` ~ ", batch_col))
    fit <- lm(form, data = df)
    out[[feat]] <- residuals(fit) + mean(df[[feat]], na.rm = TRUE)
  }
  
  return(out)
}

feature_cols <- setdiff(colnames(patient_features), c("patient_id", "source", "collection_batch"))

patient_features_adj <- residualize_by_batch(
  df = patient_features,
  feature_cols = feature_cols,
  batch_col = "collection_batch"
)

### [10] scaling and patient clustering
mat_patient <- patient_features_adj[, feature_cols, drop = FALSE]
rownames(mat_patient) <- patient_features_adj$patient_id

mat_scaled <- scale(mat_patient, center = TRUE, scale = TRUE)
mat_scaled <- as.matrix(mat_scaled)
#PCA
pca_patient <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)

pca_df <- as.data.frame(pca_patient$x[, 1:10, drop = FALSE])
pca_df$patient_id <- rownames(pca_df)
pca_df$source <- patient_features_adj$source[match(pca_df$patient_id, patient_features_adj$patient_id)]

ggplot(pca_df, aes(PC1, PC2, color = source)) +
  geom_point(size = 3) +
  theme_classic()
#UMAP
seu_patient <- CreateSeuratObject(
  counts = t(mat_scaled),
  assay = "patient"
)

seu_patient <- ScaleData(seu_patient, verbose = FALSE)
seu_patient <- RunPCA(seu_patient, npcs = 20, verbose = FALSE)
seu_patient <- RunUMAP(seu_patient, dims = 1:10)
seu_patient <- FindNeighbors(seu_patient, dims = 1:10)
seu_patient <- FindClusters(seu_patient, resolution = 0.3)
# add to metadata
seu_patient$patient_id <- colnames(seu_patient)
seu_patient$source <- patient_features_adj$source[match(colnames(seu_patient), patient_features_adj$patient_id)]
seu_patient$collection_batch <- patient_features_adj$collection_batch[match(colnames(seu_patient), patient_features_adj$patient_id)]
DimPlot(seu_patient, group.by = "source", pt.size = 2)
DimPlot(seu_patient, group.by = "seurat_clusters", label = TRUE, pt.size = 2)

### [10] hierarchical clustering
dist_mat <- dist(mat_scaled, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
plot(hc, cex = 0.7)

k <- 3
patient_cluster <- cutree(hc, k = k)
cluster_df <- data.frame(
  patient_id = names(patient_cluster),
  patient_cluster = as.factor(patient_cluster)
)
patient_features_adj <- left_join(patient_features_adj, cluster_df, by = "patient_id")

### [11] check whether normal act as negative control
table(patient_features_adj$source, patient_features_adj$patient_cluster)
DimPlot(seu_patient, group.by = "source", shape.by = "source", pt.size = 3)

