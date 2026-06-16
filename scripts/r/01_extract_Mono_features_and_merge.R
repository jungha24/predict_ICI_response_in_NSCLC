# 2026.03.25
# make feature library in Mono

# library
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(tidyr)
library(scales)
library(tidytext)

# data 
load('data/20260309_pilot/results/20260310_integrated.RData')
selected_patients <- fread('data/20260309_pilot/results/20260313_patients_selected.txt',header=F)
meta <- fread('data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt')
## meta
idx <- match(seu_integrated@meta.data$sample, meta$sample_id)
seu_integrated$binarized_response <- meta$`Binarized response`[idx]
seu_integrated$`PD.L1_TPS` <- meta$`PD-L1_TPS_mod`[idx]

# check 
DimPlot(seu_integrated, reduction = "umap", group.by='manual.cluster',split.by='manual.cluster')+facet_wrap(~manual.cluster)

# Mono subclustering

Mono_subset<- subset(seu_integrated, idents = c('CD14 Mono','CD16 Mono'))
DefaultAssay(Mono_subset) <- 'RNA'
Mono_subset$source_obj <- ifelse(grepl("^DIS_", colnames(Mono_subset)), "disease", "normal")
Mono_list <- SplitObject(Mono_subset, split.by = "source_obj")
names(Mono_list)

Mono_list <- lapply(Mono_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) #<<<
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000) #<<<
  return(x)
})

features <- SelectIntegrationFeatures(object.list = Mono_list, nfeatures = 3000) #<<<

Mono_list <- lapply(Mono_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

anchors <- FindIntegrationAnchors(
  object.list = Mono_list,
  anchor.features = features,
  reduction = "rpca", #<<<
  dims = 1:30 #<<<
)

subset <- IntegrateData(
  anchorset = anchors,
  dims = 1:30 #<<<
)

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = TRUE)
DimPlot(subset, reduction='pca', group.by='predicted.celltype.l2')
DimPlot(subset, reduction='pca', group.by='manual.cluster_l2')
ElbowPlot(subset) #3x4 #20260324_B_subcluster_PCA_elbowplot
DimHeatmap(subset, dims = 1:10, cells = 500, balanced = TRUE) #8x8 #20260324_B_subcluster_PCA_pcaheatmap
subset <- FindNeighbors(subset, dims = 1:10) #<<< 
subset <- FindClusters(subset, resolution = 0.1) #<<<
subset <- RunUMAP(subset, dims = 1:10) #<<<
DimPlot(subset, reduction='umap')#3x4
DimPlot(subset, reduction='umap', group.by = 'predicted.celltype.l2',label=T)#3x4
DimPlot(subset(subset, predicted.celltype.l2 %in% c('B intermediate','B memory','B naive')), reduction='umap',group.by='predicted.celltype.l2', split.by = 'predicted.celltype.l2')#4x8 #20260324_b_subset_cite_l2
markers <- FindAllMarkers(subset, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

subset$manual.cluster_l2 <- ifelse(subset$seurat_clusters==0,'platelet-like contam',
                                   ifelse(subset$seurat_clusters==1, 'activated monocyte',
                                          ifelse(subset$seurat_clusters==2,'CD16 monotype',
                                                 'Tcell contam')))


NONCLASSIC_LABEL = ["FCGR3A","CDKN1C","SIGLEC10","CSF1R","MTSS1","RHOC","CSF1R"]
ACTIVE_LABEL = ["IL1B","CXCL8","CCL3","TNFAIP3","IER3","LYZ"]
# CSF1R : NON-CLASSIC
# LYZ : CLASSIC
VlnPlot(subset_clean, features=c("IL1B","CXCL8","CCL3","ATF3","NR4A2","TNFAIP3","IER3","PDE4B","LYZ"),pt.size=0)
FeaturePlot(subset_clean, features=c("IL1B","CXCL8","CCL3","ATF3","NR4A2","TNFAIP3","IER3","PDE4B"), order=T)
ggplot(subset@meta.data, aes(x=manual.cluster_l2, fill=cohort_l2))+geom_bar(position='dodge')
DimPlot(subset, reduction='umap', group.by='manual.cluster_l2',label=T)#6x4 #20260326_Mono_subcluster
saveRDS(subset, file='data/20260309_pilot/results/version2/Monocyte/Mono_subset.rds')


#########
bad_clusters <- c("platelet-like contam", "Tcell contam")
subset_clean <- subset(subset, subset = !(manual.cluster_l2 %in% bad_clusters))
DefaultAssay(subset_clean) <- 'RNA'
subset_subset_list <- SplitObject(subset_clean, split.by = "source_obj")

subset_subset_list <- lapply(subset_subset_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) #<<<
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000) #<<<
  return(x)
})

features <- SelectIntegrationFeatures(object.list = subset_subset_list, nfeatures = 3000) #<<<

subset_subset_list <- lapply(subset_subset_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

anchors <- FindIntegrationAnchors(
  object.list = subset_subset_list,
  anchor.features = features,
  reduction = "rpca", #<<<
  dims = 1:30 #<<<
)

subset_clean <- IntegrateData(
  anchorset = anchors,
  dims = 1:30 #<<<
)

subset_clean <- ScaleData(subset_clean, verbose = FALSE)
subset_clean <- RunPCA(subset_clean, npcs = 50, verbose = TRUE)

ElbowPlot(subset_clean) #3x4 #20260324_B_subcluster_PCA_elbowplot
DimHeatmap(subset_clean, dims = 10:15, cells = 500, balanced = TRUE) #8x8 #20260324_B_subcluster_PCA_pcaheatmap
subset_clean <- FindNeighbors(subset_clean, dims = 1:10) #<<< 
subset_clean <- FindClusters(subset_clean, resolution = 0.2) #<<<
subset_clean <- RunUMAP(subset_clean, dims = 1:10) #<<<
DimPlot(subset_clean, reduction='umap',label=T)#4x6 #20260327_NK_exclude_IER_subcluster
markers <- FindAllMarkers(subset_clean, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

######### recalculate PCA #########
subset_clean <- subset(subset_clean, subset = cohort_l2 == 'Lee_p1_base')
DefaultAssay(subset_clean) <- 'RNA'
subset_clean <- NormalizeData(subset_clean, normalization.method = "LogNormalize", scale.factor = 10000) 
subset_clean <- FindVariableFeatures(subset_clean, selection.method = "vst", nfeatures = 3000) #<<<

subset_clean <- ScaleData(subset_clean, verbose = FALSE)
subset_clean <- RunPCA(subset_clean, npcs = 50, verbose = TRUE)

DimPlot(subset_clean, reduction='pca', group.by='predicted.celltype.l2')
DimPlot(subset_clean, reduction='pca', group.by='manual.cluster_l2')
ElbowPlot(subset_clean) #3x4 #20260324_B_subcluster_PCA_elbowplot
DimHeatmap(subset_clean, dims = 1:10, cells = 500, balanced = TRUE) #8x8 #20260328_Mono_clean_PCA
###

DefaultAssay(subset_clean) <-'RNA'
subset_clean<- JoinLayers(subset_clean)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(tibble)
library(scales)
library(edgeR)
library(msigdbr)
library(singscore)
library(ssGSEA2)
library(GeneNMF)
library(RColorBrewer)
library(Matrix)


# =========================================================
# 0. USER CONFIG
# =========================================================
TARGET_ANALYSIS <- "Monocyte"   # "Monocyte" or "NK"

# input: clean subset that already has manual.cluster_l2
# e.g.
# Monocyte: "data/20260309_pilot/results/version2/Mono/Mono_subset_clean_manual_l2.rds"
# NK:       "data/20260309_pilot/results/version2/NK/NK_subset_clean_manual_l2.rds"
#CLEAN_SUBSET_RDS <- ""

# optional: reload clinical metadata from original meta file
ADD_CLINICAL_META <- FALSE
META_PATH <- "data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt"

# common metadata columns in your Seurat object
SAMPLE_COL        <- "sample"
COHORT_COL        <- "cohort_l2"
RESPONSE_COL      <- "binarized_response"
ASSAY_USE         <- "RNA"
TARGET_COHORT     <- "Lee_p1_base"
SUBTYPE_COL       <- "manual.cluster_l2"

# optional additional cleanup after loading clean subset
# leave empty in the usual case
EXTRA_DROP_STATES <- character(0)

# patient-level thresholds
MIN_CELLS_PROP      <- 10
MIN_CELLS_PB        <- 20
MIN_NONZERO_PER_SMP <- 10
MIN_SAMPLES_DETECT  <- 2
MIN_GENESET_OVERLAP <- 3

# centroid / anchor score
# Leave as NA to skip.
# Use labels that already exist in manual.cluster_l2 of the clean subset.
ANCHOR_STATE_LOW  <- 'activated monocyte'
ANCHOR_STATE_HIGH <- 'CD16 monotype'
CENTROID_PREFIX   <- if (TARGET_ANALYSIS == "Monocyte") "mono_low_to_high" else "nk_low_to_high"

# optional merge with existing patient feature matrix
MERGE_WITH_EXISTING_FINAL_FEATURE <- FALSE
EXISTING_FINAL_FEATURE_CSV <- ""

# optional CellRank export / merge
DO_CELLRANK_EXPORT <- TRUE
KEEP_STATES_FOR_CELLRANK <- NULL  # NULL = use all states present in the clean subset
CELLRANK_PATIENT_CSV <- "data/20260309_pilot/results/version2/Monocyte/CellRank_output_from_clean_subset/cellrank_b3_patient_features.csv"       # if provided and exists, merge CellRank patient features

# output root
OUTROOT <- file.path("data/20260309_pilot/results/version2", TARGET_ANALYSIS)
dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)
# =========================================================
# 1. ANALYSIS-SPECIFIC DEFAULTS
# =========================================================
get_analysis_defaults <- function(target_analysis) {
  if (target_analysis == "Monocyte") {
    list(
      gene_blocklist_patterns = c("^MT-", "^RPS", "^RPL", "^HSP", "^IG[HKL]", "^TRA", "^TRB", "^TRD", "^TRG"),
      custom_gene_sets = list(
        MONO_MHCII_AP = c("CD74", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "CIITA", "IFI30", "LGMN", "CTSB", "CTSD", "CTSL", "CTSS"),
        MONO_IFN_CORE = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI27", "IFI44", "IFI44L", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "STAT1", "IRF7", "IFITM3", "CXCL10"),
        MONO_INFLAMMATORY = c("IL1B", "TNF", "NFKBIA", "CXCL8", "CCL3", "CCL4", "S100A8", "S100A9", "PTGS2", "BCL2A1", "EREG", "FOS", "JUNB"),
        MONO_FCGR3A_MYELOID = c("FCGR3A", "LST1", "IFITM3", "SAT1", "TYMP", "CTD", "AIF1", "MS4A7", "COTL1", "CFD", "RHOC", "LGALS3")
      ),
      hallmark_targets = c(
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
        "GOBP_MONOCYTE_CHEMOTAXIS",
        "GOBP_REGULATION_OF_MONOCYTE_CHEMOTAXIS",
        "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
        "GOBP_PHAGOCYTOSIS",
        "GOBP_PHAGOCYTOSIS_RECOGNITION",
        "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
        "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE",
        "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
        "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
        "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT"
      )
    )
  } else if (target_analysis == "NK") {
    list(
      gene_blocklist_patterns = c("^MT-", "^RPS", "^RPL", "^HSP", "^IG[HKL]", "^TRA", "^TRB", "^TRD", "^TRG"),
      custom_gene_sets = list(
        NK_CYTOTOXICITY = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "CTSW", "CST7", "FGFBP2", "KLRD1", "FCGR3A", "CCL5", "XCL1", "XCL2"),
        NK_IFNG_RESPONSE = c("IFNG", "STAT1", "IRF1", "CXCL10", "CXCL9", "GBP1", "ISG15", "IFIT3", "OAS1", "OASL", "MX1"),
        NK_ACTIVATION = c("CD69", "TNF", "XCL1", "XCL2", "IFNG", "NFKBIA", "JUN", "FOS", "DUSP1", "NR4A1", "NR4A2"),
        NK_ADAPTIVE_FCGR3A = c("FCGR3A", "TYROBP", "KLRC2", "LILRB1", "ZEB2", "B3GAT1", "SPON2", "TRAC", "KLRD1", "CCL5")
      ),
      hallmark_targets = c(
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
        "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
        "GOBP_NATURAL_KILLER_CELL_CHEMOTAXIS",
        "GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
        "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
        "GOBP_POSITIVE_REGULATION_OF_TYPE_II_INTERFERON_PRODUCTION"
      )
    )
  } else {
    stop("TARGET_ANALYSIS must be 'Monocyte' or 'NK'.")
  }
}

cfg <- get_analysis_defaults(TARGET_ANALYSIS)
GENE_BLOCKLIST_PATTERNS <- cfg$gene_blocklist_patterns
CUSTOM_GENE_SETS <- cfg$custom_gene_sets
# HALLMARK_TARGETS <- cfg$hallmark_targets
MSIGDB_TARGETS <- cfg$hallmark_targets
# =========================================================
# 2. HELPERS
# =========================================================
`%||%` <- function(x, y) if (!is.null(x)) x else y

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

standardize_sample_for_merge <- function(x) {
  x <- as.character(x)
  gsub("_", "-", x)
}

check_clean_subset <- function(seu, subtype_col, sample_col, cohort_col, assay_use) {
  stopifnot(inherits(seu, "Seurat"))
  if (!subtype_col %in% colnames(seu@meta.data)) {
    stop(paste0("Column not found in metadata: ", subtype_col))
  }
  if (!sample_col %in% colnames(seu@meta.data)) {
    stop(paste0("Column not found in metadata: ", sample_col))
  }
  if (!cohort_col %in% colnames(seu@meta.data)) {
    stop(paste0("Column not found in metadata: ", cohort_col))
  }
  if (!assay_use %in% Assays(seu)) {
    stop(paste0("Assay not found: ", assay_use))
  }
  if (all(is.na(seu@meta.data[[subtype_col]]))) {
    stop(paste0(subtype_col, " exists but is all NA."))
  }
}

calc_two_anchor_score <- function(
    seu,
    sample_col = "sample",
    celltype_col = "manual.cluster_l2",
    low_label,
    high_label,
    assay_use = "RNA",
    rerun_pca = TRUE,
    reduction_use = "pca",
    dims_use = c(1,4,5,7),
    npcs = 20,
    min_cells_per_patient = 10,
    score_cut = 0.5,
    eps = 1e-8,
    prefix = "state_score"
) {
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(sample_col %in% colnames(seu@meta.data))
  stopifnot(celltype_col %in% colnames(seu@meta.data))
  
  meta_df <- seu@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    mutate(
      .sample = .data[[sample_col]],
      .celltype = .data[[celltype_col]]
    )
  
  if (!all(c(low_label, high_label) %in% unique(meta_df$.celltype))) {
    stop("low_label/high_label not found in celltype_col.")
  }
  
  if (sum(meta_df$.celltype == low_label) < 2) stop("Too few low anchor cells.")
  if (sum(meta_df$.celltype == high_label) < 2) stop("Too few high anchor cells.")
  
  seu_work <- seu
  DefaultAssay(seu_work) <- assay_use
  
  if (rerun_pca) {
    seu_work <- NormalizeData(seu_work, verbose = FALSE)
    seu_work <- FindVariableFeatures(seu_work, verbose = FALSE)
    seu_work <- ScaleData(seu_work, verbose = FALSE)
    seu_work <- RunPCA(seu_work, npcs = npcs, verbose = FALSE)
    reduction_use <- "pca"
  } else {
    if (!reduction_use %in% Reductions(seu_work)) {
      stop(paste0("Reduction '", reduction_use, "' not found."))
    }
  }
  
  emb <- Embeddings(seu_work, reduction = reduction_use)
  if (max(dims_use) > ncol(emb)) stop("dims_use exceeds available dimensions.")
  emb <- emb[, dims_use, drop = FALSE]
  
  df_emb <- emb %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    inner_join(meta_df, by = "cell")
  
  dim_cols <- colnames(emb)
  low_centroid <- colMeans(df_emb[df_emb$.celltype == low_label, dim_cols, drop = FALSE])
  high_centroid <- colMeans(df_emb[df_emb$.celltype == high_label, dim_cols, drop = FALSE])
  
  v <- high_centroid - low_centroid
  vv <- sum(v^2)
  if (vv == 0) stop("Low and high centroids are identical.")
  
  X <- as.matrix(df_emb[, dim_cols, drop = FALSE])
  raw_score <- as.numeric(
    (X - matrix(low_centroid, nrow = nrow(X), ncol = length(low_centroid), byrow = TRUE)) %*% v / vv
  )
  clipped_score <- pmin(pmax(raw_score, 0), 1)
  
  score_col     <- paste0(prefix, "_score")
  raw_col       <- paste0(prefix, "_score_raw")
  bin_col       <- paste0(prefix, "_bin")
  low_frac_col  <- paste0("frac_", sanitize_name(low_label), "_like")
  high_frac_col <- paste0("frac_", sanitize_name(high_label), "_like")
  obs_low_col   <- paste0("observed_frac_", sanitize_name(low_label))
  obs_high_col  <- paste0("observed_frac_", sanitize_name(high_label))
  logratio_col  <- paste0("log2ratio_", sanitize_name(high_label), "_", sanitize_name(low_label))
  skew_col      <- paste0("skew_", sanitize_name(high_label), "_minus_", sanitize_name(low_label))
  
  df_emb[[raw_col]] <- raw_score
  df_emb[[score_col]] <- clipped_score
  df_emb[[bin_col]] <- ifelse(df_emb[[score_col]] < score_cut, paste0(low_label, "_like"), paste0(high_label, "_like"))
  
  summary_df <- df_emb %>%
    group_by(.sample) %>%
    summarise(
      n_score_cells = n(),
      mean_score = mean(.data[[score_col]], na.rm = TRUE),
      median_score = median(.data[[score_col]], na.rm = TRUE),
      sd_score = sd(.data[[score_col]], na.rm = TRUE),
      iqr_score = IQR(.data[[score_col]], na.rm = TRUE),
      !!low_frac_col  := mean(.data[[bin_col]] == paste0(low_label, "_like"), na.rm = TRUE),
      !!high_frac_col := mean(.data[[bin_col]] == paste0(high_label, "_like"), na.rm = TRUE),
      !!obs_low_col   := mean(.celltype == low_label, na.rm = TRUE),
      !!obs_high_col  := mean(.celltype == high_label, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      !!logratio_col := log2((.data[[high_frac_col]] + eps) / (.data[[low_frac_col]] + eps)),
      !!skew_col := .data[[high_frac_col]] - .data[[low_frac_col]],
      keep_for_analysis = n_score_cells >= min_cells_per_patient
    ) %>%
    arrange(median_score)
  
  meta_add <- df_emb %>%
    select(cell, all_of(c(score_col, raw_col, bin_col)))
  meta_add2 <- meta_add[match(Cells(seu), meta_add$cell), -1, drop = FALSE]
  rownames(meta_add2) <- Cells(seu)
  seu_scored <- AddMetaData(seu, metadata = meta_add2)
  
  list(
    seu_scored = seu_scored,
    cell_scores = df_emb,
    patient_summary = summary_df,
    low_centroid = low_centroid,
    high_centroid = high_centroid,
    prefix = prefix,
    dims_used = dims_use,
    reduction_used = reduction_use
  )
}

make_patient_composition <- function(seu, sample_col, cohort_col, response_col, celltype_col) {
  meta_df <- seu@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  
  meta_df[[sample_col]] <- ifelse(
    is.na(meta_df[[sample_col]]),
    meta_df$donor_id %||% meta_df[[sample_col]],
    meta_df[[sample_col]]
  )
  
  if (response_col %in% colnames(meta_df)) {
    meta_df$cohort_l2_resp <- paste0(meta_df[[cohort_col]], "_", meta_df[[response_col]])
  } else {
    meta_df$cohort_l2_resp <- meta_df[[cohort_col]]
  }
  
  comp_long <- meta_df %>%
    filter(!is.na(.data[[sample_col]]), .data[[sample_col]] != "") %>%
    count(cohort_l2_resp, .data[[sample_col]], .data[[celltype_col]], name = "n_cells") %>%
    rename(sample = all_of(sample_col), manual.cluster_l2 = all_of(celltype_col)) %>%
    group_by(cohort_l2_resp, sample) %>%
    mutate(total_cells = sum(n_cells), prop = n_cells / total_cells) %>%
    ungroup()
  
  comp_long
}

run_curated_pseudobulk_scores <- function(
    seu_use,
    outdir,
    sample_col,
    cohort_col = NULL,
    target_cohort = NULL,
    assay_use = "RNA",
    species_use = "Homo sapiens",
    gene_sets = list(),
    hallmark_targets = NULL,
    min_cells_per_sample = 10,
    min_nonzero_cells_per_sample = 1,
    min_samples_passing_detection = 2,
    min_geneset_overlap = 3,
    prefix_name = "curated"
) {
  ## ---------------------------------------------------------
  ## packages / checks
  ## ---------------------------------------------------------
  stopifnot(inherits(seu_use, "Seurat"))
  
  req_pkgs <- c("SeuratObject", "Matrix", "dplyr", "tibble", "msigdbr")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  meta <- seu_use@meta.data
  if (!(sample_col %in% colnames(meta))) {
    stop("sample_col not found in Seurat metadata: ", sample_col)
  }
  if (!(assay_use %in% names(seu_use@assays))) {
    stop("assay_use not found in Seurat object: ", assay_use)
  }
  
  ## ---------------------------------------------------------
  ## helpers
  ## ---------------------------------------------------------
  species_lower <- tolower(species_use)
  is_mouse <- species_lower %in% c("mouse", "mus musculus")
  
  norm_symbol <- function(x) {
    x <- as.character(x)
    if (is_mouse) x else toupper(x)
  }
  
  empty_return <- function(overlap_tbl = tibble::tibble(), gene_sets_filt = list()) {
    empty_ssgsea <- tibble::tibble(
      sample = character(),
      gs_name = character(),
      ssgsea2_es = numeric()
    )
    empty_singscore <- tibble::tibble(
      sample = character(),
      gs_name = character(),
      singscore = numeric()
    )
    empty_pc <- tibble::tibble(
      sample = character(),
      gs_name = character(),
      pc1 = numeric(),
      pc2 = numeric(),
      eigengene = numeric()
    )
    empty_feature <- tibble::tibble(
      sample = character(),
      gs_name = character()
    )
    
    return(list(
      overlap_tbl = overlap_tbl,
      gene_sets_filt = gene_sets_filt,
      pb_counts = NULL,
      pb_cpm = NULL,
      pb_logcpm = NULL,
      ssgsea_long = empty_ssgsea,
      singscore_long = empty_singscore,
      pc_long = empty_pc,
      feature_partial = empty_feature,
      feature_long = empty_feature
    ))
  }
  
  ## ---------------------------------------------------------
  ## 1. cell filtering
  ## ---------------------------------------------------------
  cells_use <- colnames(seu_use)
  
  if (!is.null(cohort_col) && !is.null(target_cohort)) {
    if (!(cohort_col %in% colnames(meta))) {
      stop("cohort_col not found in Seurat metadata: ", cohort_col)
    }
    keep_idx <- meta[[cohort_col]] %in% target_cohort
    cells_use <- rownames(meta)[keep_idx]
  }
  
  if (length(cells_use) == 0) {
    warning("No cells left after cohort filtering.")
    return(empty_return())
  }
  
  meta_use <- meta[cells_use, , drop = FALSE]
  sample_vec <- as.character(meta_use[[sample_col]])
  names(sample_vec) <- rownames(meta_use)
  
  sample_ncells <- sort(table(sample_vec), decreasing = TRUE)
  keep_samples <- names(sample_ncells)[sample_ncells >= min_cells_per_sample]
  
  message("\n[INFO] samples before min_cells filter: ", length(sample_ncells))
  message("[INFO] samples after min_cells filter: ", length(keep_samples))
  
  if (length(keep_samples) == 0) {
    warning("No samples passed min_cells_per_sample = ", min_cells_per_sample)
    return(empty_return())
  }
  
  cells_keep <- names(sample_vec)[sample_vec %in% keep_samples]
  meta_keep <- meta_use[cells_keep, , drop = FALSE]
  sample_vec_keep <- as.character(meta_keep[[sample_col]])
  names(sample_vec_keep) <- rownames(meta_keep)
  
  ## ---------------------------------------------------------
  ## 2. counts extraction
  ## ---------------------------------------------------------
  counts_use <- LayerData(
    seu_use,
    assay = assay_use,
    layer = "counts"
  )[, cells_keep, drop = FALSE]
  
  if (ncol(counts_use) == 0) {
    warning("No cells remain after filtering.")
    return(empty_return())
  }
  
  ## ---------------------------------------------------------
  ## 3. pseudobulk counts (gene x sample)
  ## ---------------------------------------------------------
  sample_fac <- factor(sample_vec_keep, levels = keep_samples)
  sample_mm <- Matrix::sparse.model.matrix(~ 0 + sample_fac)
  colnames(sample_mm) <- levels(sample_fac)
  
  pb_counts <- counts_use %*% sample_mm
  pb_counts <- as.matrix(pb_counts)
  
  message("\n[INFO] dim(pb_counts): ", paste(dim(pb_counts), collapse = " x "))
  
  if (ncol(pb_counts) == 0) {
    warning("pb_counts has 0 columns.")
    return(empty_return())
  }
  
  ## ---------------------------------------------------------
  ## 4. CPM / logCPM
  ## ---------------------------------------------------------
  lib_sizes <- colSums(pb_counts)
  if (any(lib_sizes == 0)) {
    zero_lib <- names(lib_sizes)[lib_sizes == 0]
    warning("Dropping samples with zero library size: ", paste(zero_lib, collapse = ", "))
    keep_lib <- lib_sizes > 0
    pb_counts <- pb_counts[, keep_lib, drop = FALSE]
    lib_sizes <- lib_sizes[keep_lib]
  }
  
  if (ncol(pb_counts) == 0) {
    warning("No samples remain after removing zero-library pseudobulk samples.")
    return(empty_return())
  }
  
  pb_cpm <- t(t(pb_counts) / lib_sizes) * 1e6
  pb_logcpm <- log2(pb_cpm + 1)
  
  ## ---------------------------------------------------------
  ## 5. gene-level nonzero cell detection table (gene x sample)
  ## ---------------------------------------------------------
  nonzero_by_gene_sample <- as.matrix((counts_use > 0) %*% sample_mm)
  
  rownames(nonzero_by_gene_sample) <- norm_symbol(rownames(nonzero_by_gene_sample))
  
  message("[INFO] dim(nonzero_by_gene_sample): ",
          paste(dim(nonzero_by_gene_sample), collapse = " x "))
  
  ## ---------------------------------------------------------
  ## 6. prepare gene sets
  ##    - supports mixed targets:
  ##      HALLMARK_*
  ##      GOBP_*
  ##      REACTOME_*
  ##      KEGG_*
  ## ---------------------------------------------------------
  load_msigdb_targets <- function(targets, species_use, norm_symbol) {
    targets <- unique(as.character(targets))
    targets <- targets[!is.na(targets) & targets != ""]
    if (length(targets) == 0) return(list())
    
    targets_h  <- targets[grepl("^HALLMARK_", targets)]
    targets_bp <- targets[grepl("^GOBP_", targets)]
    targets_re <- targets[grepl("^REACTOME_", targets)]
    targets_kg <- targets[grepl("^KEGG_", targets)]
    
    out <- list()
    
    ## Hallmark
    if (length(targets_h) > 0) {
      df_h <- msigdbr::msigdbr(
        species = species_use,
        collection = "H"
      ) %>%
        dplyr::filter(.data$gs_name %in% targets_h) %>%
        dplyr::transmute(
          gs_name = .data$gs_name,
          gene_symbol = norm_symbol(.data$gene_symbol)
        ) %>%
        dplyr::distinct()
      
      sets_h <- split(df_h$gene_symbol, df_h$gs_name)
      sets_h <- sets_h[intersect(targets_h, names(sets_h))]
      out <- c(out, sets_h)
    }
    
    ## GO biological process
    if (length(targets_bp) > 0) {
      df_bp <- msigdbr::msigdbr(
        species = species_use,
        collection = "C5",
        subcollection = "GO:BP"
      ) %>%
        dplyr::filter(.data$gs_name %in% targets_bp) %>%
        dplyr::transmute(
          gs_name = .data$gs_name,
          gene_symbol = norm_symbol(.data$gene_symbol)
        ) %>%
        dplyr::distinct()
      
      sets_bp <- split(df_bp$gene_symbol, df_bp$gs_name)
      sets_bp <- sets_bp[intersect(targets_bp, names(sets_bp))]
      out <- c(out, sets_bp)
    }
    
    ## Reactome
    if (length(targets_re) > 0) {
      df_re <- msigdbr::msigdbr(
        species = species_use,
        collection = "C2",
        subcollection = "CP:REACTOME"
      ) %>%
        dplyr::filter(.data$gs_name %in% targets_re) %>%
        dplyr::transmute(
          gs_name = .data$gs_name,
          gene_symbol = norm_symbol(.data$gene_symbol)
        ) %>%
        dplyr::distinct()
      
      sets_re <- split(df_re$gene_symbol, df_re$gs_name)
      sets_re <- sets_re[intersect(targets_re, names(sets_re))]
      out <- c(out, sets_re)
    }
    
    ## KEGG
    ## try both legacy + medicus, because current msigdbr separates them
    if (length(targets_kg) > 0) {
      df_kg_legacy <- msigdbr::msigdbr(
        species = species_use,
        collection = "C2",
        subcollection = "CP:KEGG_LEGACY"
      ) %>%
        dplyr::filter(.data$gs_name %in% targets_kg) %>%
        dplyr::transmute(
          gs_name = .data$gs_name,
          gene_symbol = norm_symbol(.data$gene_symbol)
        ) %>%
        dplyr::distinct()
      
      df_kg_medicus <- msigdbr::msigdbr(
        species = species_use,
        collection = "C2",
        subcollection = "CP:KEGG_MEDICUS"
      ) %>%
        dplyr::filter(.data$gs_name %in% targets_kg) %>%
        dplyr::transmute(
          gs_name = .data$gs_name,
          gene_symbol = norm_symbol(.data$gene_symbol)
        ) %>%
        dplyr::distinct()
      
      df_kg <- dplyr::bind_rows(df_kg_legacy, df_kg_medicus) %>%
        dplyr::distinct()
      
      sets_kg <- split(df_kg$gene_symbol, df_kg$gs_name)
      sets_kg <- sets_kg[intersect(targets_kg, names(sets_kg))]
      out <- c(out, sets_kg)
    }
    
    out
  }
  
  msigdb_sets <- list()
  
  if (!is.null(hallmark_targets) && length(hallmark_targets) > 0) {
    msigdb_sets <- load_msigdb_targets(
      targets = hallmark_targets,
      species_use = species_use,
      norm_symbol = norm_symbol
    )
  }
  
  custom_sets <- gene_sets
  if (length(custom_sets) > 0) {
    custom_sets <- lapply(custom_sets, function(x) unique(norm_symbol(x)))
  }
  
  gene_sets_all <- c(msigdb_sets, custom_sets)
  
  if (length(gene_sets_all) == 0) {
    warning("No MSigDB targets or custom gene_sets supplied.")
    return(empty_return())
  }
  
  message("\n[INFO] total input gene sets:")
  print(names(gene_sets_all))
  
  missing_targets <- setdiff(
    unique(as.character(hallmark_targets)),
    names(msigdb_sets)
  )
  if (length(missing_targets) > 0) {
    message("\n[INFO] MSigDB targets not found in selected collections/subcollections:")
    print(missing_targets)
  }
  
  ## ---------------------------------------------------------
  ## 7. gene set-level detection filtering
  ##    - do NOT globally filter pb_logcpm
  ##    - filter only within each gene set
  ## ---------------------------------------------------------
  expr_genes <- rownames(pb_logcpm)
  expr_genes_norm <- norm_symbol(expr_genes)
  
  filter_one_geneset <- function(gs) {
    gs <- unique(norm_symbol(gs))
    
    ## genes available in detection table
    gs0 <- intersect(gs, rownames(nonzero_by_gene_sample))
    if (length(gs0) == 0) return(character(0))
    
    det_mat <- nonzero_by_gene_sample[gs0, , drop = FALSE]
    n_samples_pass <- rowSums(det_mat >= min_nonzero_cells_per_sample)
    gs1 <- gs0[n_samples_pass >= min_samples_passing_detection]
    
    ## map back to original expr gene names
    expr_genes[expr_genes_norm %in% gs1]
  }
  
  gene_sets_detect_filt <- lapply(gene_sets_all, filter_one_geneset)
  
  overlap_tbl <- tibble::tibble(
    gs_name = names(gene_sets_all),
    n_input_genes = vapply(gene_sets_all, length, integer(1)),
    n_after_detection = vapply(gene_sets_detect_filt, length, integer(1))
  ) %>%
    dplyr::mutate(
      n_overlap = .data$n_after_detection,
      status = dplyr::case_when(
        .data$n_after_detection == 0 ~ "skip_detection",
        .data$n_after_detection < min_geneset_overlap ~ "skip_overlap",
        TRUE ~ "keep"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_overlap), .data$gs_name)
  
  write.csv(
    overlap_tbl,
    file.path(outdir, paste0(prefix_name, "_06_gene_set_overlap.csv")),
    row.names = FALSE
  )
  
  message("\n[INFO] gene set overlap / status:")
  print(overlap_tbl)
  
  keep_gs <- overlap_tbl %>%
    dplyr::filter(.data$status == "keep") %>%
    dplyr::pull(.data$gs_name)
  
  gene_sets_filt <- gene_sets_detect_filt[keep_gs]
  
  saveRDS(
    gene_sets_filt,
    file.path(outdir, paste0(prefix_name, "_06_gene_sets_filtered.rds"))
  )
  
  message("\n[INFO] final gene sets kept:")
  print(names(gene_sets_filt))
  
  if (length(gene_sets_filt) == 0) {
    warning(
      "No gene sets passed detection/overlap filtering. ",
      "Returning empty result."
    )
    return(empty_return(overlap_tbl = overlap_tbl, gene_sets_filt = gene_sets_filt))
  }
  
  ## ---------------------------------------------------------
  ## 8. export genes = union of kept gene sets
  ## ---------------------------------------------------------
  genes_for_export <- unique(unlist(gene_sets_filt, use.names = FALSE))
  pb_logcpm_export <- pb_logcpm[genes_for_export, , drop = FALSE]
  
  message("\n[INFO] dim(pb_logcpm_export): ",
          paste(dim(pb_logcpm_export), collapse = " x "))
  
  if (nrow(pb_logcpm_export) == 0 || ncol(pb_logcpm_export) == 0) {
    warning("pb_logcpm_export is empty after gene set filtering.")
    return(empty_return(overlap_tbl = overlap_tbl, gene_sets_filt = gene_sets_filt))
  }
  
  ## ---------------------------------------------------------
  ## 9. GMT export
  ## ---------------------------------------------------------
  gmt_file <- file.path(outdir, paste0(prefix_name, "_gene_sets.gmt"))
  
  gmt_lines <- vapply(names(gene_sets_filt), function(gs) {
    paste(c(gs, "na", gene_sets_filt[[gs]]), collapse = "\t")
  }, character(1))
  
  writeLines(gmt_lines, con = gmt_file)
  
  ## ---------------------------------------------------------
  ## 10. GCT export
  ## ---------------------------------------------------------
  gct_file <- file.path(outdir, paste0(prefix_name, "_expr.gct"))
  
  con <- file(gct_file, open = "wt")
  writeLines("#1.2", con)
  writeLines(sprintf("%d\t%d", nrow(pb_logcpm_export), ncol(pb_logcpm_export)), con)
  writeLines(
    paste(c("Name", "Description", colnames(pb_logcpm_export)), collapse = "\t"),
    con
  )
  
  gct_df <- data.frame(
    Name = rownames(pb_logcpm_export),
    Description = rep("na", nrow(pb_logcpm_export)),
    pb_logcpm_export,
    check.names = FALSE
  )
  
  write.table(
    gct_df,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  close(con)
  
  ## ---------------------------------------------------------
  ## 11. ssGSEA2
  ## ---------------------------------------------------------
  ssgsea_long <- tibble::tibble(
    sample = character(),
    gs_name = character(),
    ssgsea2_es = numeric()
  )
  
  ssgsea_raw <- NULL
  ssgsea_named <- NULL
  
  if (requireNamespace("ssGSEA2", quietly = TRUE)) {
    ssgsea_out_prefix <- paste0(prefix_name, "_ssgsea2")
    ssgsea_log_file <- file.path(outdir, paste0(ssgsea_out_prefix, ".run.log"))
    
    ssgsea_ok <- TRUE
    ssgsea_err <- NULL
    
    tryCatch({
      ssgsea_raw <- ssGSEA2::run_ssGSEA2(
        input.ds = gct_file,
        output.prefix = ssgsea_out_prefix,
        gene.set.databases = gmt_file,
        output.directory = outdir,
        sample.norm.type = "none",
        weight = 0.75,
        correl.type = "rank",
        statistic = "area.under.RES",
        output.score.type = "NES",
        nperm = 1000,
        min.overlap = min_geneset_overlap,
        extended.output = TRUE,
        global.fdr = FALSE,
        export.signat.gct = TRUE,
        param.file = TRUE,
        log.file = ssgsea_log_file
      )
    }, error = function(e) {
      ssgsea_ok <<- FALSE
      ssgsea_err <<- conditionMessage(e)
    })
    
    if (ssgsea_ok && !is.null(ssgsea_raw)) {
      sample_names <- colnames(pb_logcpm_export)
      
      ssgsea_named <- lapply(ssgsea_raw, function(pathway_obj) {
        n_use <- min(length(pathway_obj), length(sample_names))
        names(pathway_obj)[seq_len(n_use)] <- sample_names[seq_len(n_use)]
        pathway_obj
      })
      
      saveRDS(
        ssgsea_named,
        file.path(outdir, paste0(prefix_name, "_07_ssgsea2_named_result.rds"))
      )
      
      ssgsea_long <- dplyr::bind_rows(lapply(names(ssgsea_named), function(gs_name) {
        pathway_obj <- ssgsea_named[[gs_name]]
        
        tibble::tibble(
          sample = names(pathway_obj),
          gs_name = gs_name,
          ssgsea2_es = vapply(pathway_obj, function(x) {
            if (!is.null(x$ES)) as.numeric(x$ES[1]) else NA_real_
          }, numeric(1))
        )
      }))
      
      write.csv(
        ssgsea_long,
        file.path(outdir, paste0(prefix_name, "_07_ssgsea2_long.csv")),
        row.names = FALSE
      )
    } else {
      warning("ssGSEA2 failed and will be skipped. Message: ", ssgsea_err)
    }
  } else {
    warning("Package 'ssGSEA2' is not installed. Skipping ssGSEA2.")
  }
  
  ## ---------------------------------------------------------
  ## 12. singscore
  ## ---------------------------------------------------------
  singscore_long <- tibble::tibble(
    sample = character(),
    gs_name = character(),
    singscore = numeric()
  )
  
  if (requireNamespace("singscore", quietly = TRUE)) {
    rank_data <- singscore::rankGenes(pb_logcpm, tiesMethod = "min")
    
    singscore_long_list <- lapply(names(gene_sets_filt), function(gs_name) {
      ss <- singscore::simpleScore(
        rankData = rank_data,
        upSet = gene_sets_filt[[gs_name]],
        knownDirection = TRUE
      )
      
      score_col <- if ("TotalScore" %in% colnames(ss)) {
        "TotalScore"
      } else if ("UpScore" %in% colnames(ss)) {
        "UpScore"
      } else {
        colnames(ss)[1]
      }
      
      tibble::tibble(
        sample = rownames(ss),
        gs_name = gs_name,
        singscore = ss[[score_col]]
      )
    })
    
    singscore_long <- dplyr::bind_rows(singscore_long_list)
    
    write.csv(
      singscore_long,
      file.path(outdir, paste0(prefix_name, "_08_singscore_long.csv")),
      row.names = FALSE
    )
  } else {
    warning("Package 'singscore' is not installed. Skipping singscore.")
  }
  
  ## ---------------------------------------------------------
  ## 13. PC1 / PC2 / eigengene
  ## ---------------------------------------------------------
  pc_long_list <- lapply(names(gene_sets_filt), function(gs_name) {
    genes <- intersect(gene_sets_filt[[gs_name]], rownames(pb_logcpm))
    X <- t(pb_logcpm[genes, , drop = FALSE])  ## sample x genes
    
    if (nrow(X) == 0) {
      return(tibble::tibble(
        sample = character(0),
        gs_name = character(0),
        pc1 = numeric(0),
        pc2 = numeric(0),
        eigengene = numeric(0)
      ))
    }
    
    if (ncol(X) == 0) {
      return(tibble::tibble(
        sample = rownames(X),
        gs_name = gs_name,
        pc1 = NA_real_,
        pc2 = NA_real_,
        eigengene = NA_real_
      ))
    }
    
    if (ncol(X) == 1) {
      pc1 <- as.numeric(scale(X[, 1]))
      pc2 <- rep(NA_real_, length(pc1))
      names(pc1) <- rownames(X)
    } else {
      pca <- prcomp(X, center = TRUE, scale. = TRUE)
      pc1 <- pca$x[, 1]
      pc2 <- if (ncol(pca$x) >= 2) pca$x[, 2] else rep(NA_real_, nrow(X))
      
      avg_expr <- rowMeans(X, na.rm = TRUE)
      cc <- suppressWarnings(stats::cor(pc1, avg_expr, use = "pairwise.complete.obs"))
      if (!is.na(cc) && cc < 0) {
        pc1 <- -pc1
      }
    }
    
    tibble::tibble(
      sample = rownames(X),
      gs_name = gs_name,
      pc1 = as.numeric(pc1),
      pc2 = as.numeric(pc2),
      eigengene = as.numeric(pc1)
    )
  })
  
  pc_long <- dplyr::bind_rows(pc_long_list)
  
  write.csv(
    pc_long,
    file.path(outdir, paste0(prefix_name, "_09_pc_long.csv")),
    row.names = FALSE
  )
  
  ## ---------------------------------------------------------
  ## 14. merge outputs
  ## ---------------------------------------------------------
  feature_partial <- singscore_long %>%
    dplyr::full_join(pc_long, by = c("sample", "gs_name")) %>%
    dplyr::left_join(overlap_tbl, by = "gs_name")
  
  feature_long <- ssgsea_long %>%
    dplyr::full_join(singscore_long, by = c("sample", "gs_name")) %>%
    dplyr::full_join(pc_long, by = c("sample", "gs_name")) %>%
    dplyr::left_join(overlap_tbl, by = "gs_name")
  
  write.csv(
    feature_partial,
    file.path(outdir, paste0(prefix_name, "_10_feature_partial_no_ssgsea.csv")),
    row.names = FALSE
  )
  
  write.csv(
    feature_long,
    file.path(outdir, paste0(prefix_name, "_10_feature_partial_yes_ssgsea.csv")),
    row.names = FALSE
  )
  
  ## ---------------------------------------------------------
  ## 15. return
  ## ---------------------------------------------------------
  return(list(
    overlap_tbl = overlap_tbl,
    gene_sets_filt = gene_sets_filt,
    pb_counts = pb_counts,
    pb_cpm = pb_cpm,
    pb_logcpm = pb_logcpm,
    pb_logcpm_export = pb_logcpm_export,
    gmt_file = gmt_file,
    gct_file = gct_file,
    ssgsea_raw = ssgsea_raw,
    ssgsea_named = ssgsea_named,
    ssgsea_long = ssgsea_long,
    singscore_long = singscore_long,
    pc_long = pc_long,
    feature_partial = feature_partial,
    feature_long = feature_long
  ))
}

run_gene_nmf_features <- function(
    seu_use,
    outdir,
    sample_col = "sample",
    cohort_col = "cohort_l2",
    target_cohort = "Lee_p1_base",
    assay_use = "RNA",
    min_cells_per_sample = 20,
    k_vec = 4:7,
    nfeatures_hvg = 2000,
    nmf_seed = 123,
    gene_blocklist_patterns = c("^MT-",
                                "^RPS", "^RPL",
                                "^HSP",
                                "^IG[HKL]",
                                "^TRA", "^TRB", "^TRD", "^TRG")
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  meta0 <- seu_use@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  
  keep_cells_cohort <- meta0$cell[!is.na(meta0[[cohort_col]]) & meta0[[cohort_col]] == target_cohort]
  seu_cohort <- subset(seu_use, cells = keep_cells_cohort)
  
  meta1 <- seu_cohort@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  keep_cells_non_na <- meta1$cell[!is.na(meta1[[sample_col]]) & meta1[[sample_col]] != ""]
  seu_non_na <- subset(seu_cohort, cells = keep_cells_non_na)
  
  meta2 <- seu_non_na@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  sample_meta2 <- meta2 %>%
    transmute(sample = as.character(.data[[sample_col]])) %>%
    count(sample, name = "n_cells") %>%
    arrange(desc(n_cells))
  keep_samples <- sample_meta2 %>% filter(n_cells >= min_cells_per_sample) %>% pull(sample)
  
  if (length(keep_samples) < 2) {
    stop("Need at least 2 samples for robust multi-sample GeneNMF analysis.")
  }
  
  keep_cells <- meta2$cell[meta2[[sample_col]] %in% keep_samples]
  seu_nm <- subset(seu_non_na, cells = keep_cells)
  obj_list <- SplitObject(seu_nm, split.by = sample_col)
  sample_sizes <- sapply(obj_list, ncol)
  obj_list <- obj_list[sample_sizes >= min_cells_per_sample]
  
  all_genes <- rownames(obj_list[[1]])
  genes_blocklist <- unique(unlist(lapply(gene_blocklist_patterns, function(p) grep(p, all_genes, value = TRUE, ignore.case = FALSE))))
  write.csv(data.frame(gene = genes_blocklist), file.path(outdir, "01_genes_blocklist.csv"), row.names = FALSE)
  
  obj_list <- lapply(obj_list, function(x) {
    DefaultAssay(x) <- assay_use
    x <- JoinLayers(x, assay = assay_use)
    if (!"data" %in% Layers(x[[assay_use]])) x <- NormalizeData(x, assay = assay_use, verbose = FALSE)
    x
  })
  
  obj_list_hvg <- lapply(obj_list, function(x) {
    DefaultAssay(x) <- assay_use
    GeneNMF::findVariableFeatures_wfilters(
      obj = x,
      nfeatures = nfeatures_hvg,
      genesBlockList = genes_blocklist,
      min.exp = 0.01,
      max.exp = 3
    )
  })
  
  hvg_list <- lapply(obj_list_hvg, VariableFeatures)
  hvg_summary <- data.frame(sample = names(hvg_list), n_hvg = sapply(hvg_list, length))
  write.csv(hvg_summary, file.path(outdir, "02_hvg_summary.csv"), row.names = FALSE)
  
  hvg_union <- sort(unique(unlist(hvg_list)))
  write.csv(data.frame(gene = hvg_union), file.path(outdir, "02_hvg_union.csv"), row.names = FALSE)
  
  geneNMF_programs <- GeneNMF::multiNMF(
    obj.list = obj_list_hvg,
    assay = assay_use,
    slot = "data",
    k = k_vec,
    hvg = hvg_union,
    nfeatures = nfeatures_hvg,
    L1 = c(0, 0),
    min.exp = 0.01,
    max.exp = 3,
    center = FALSE,
    scale = FALSE,
    min.cells.per.sample = min_cells_per_sample,
    hvg.blocklist = genes_blocklist,
    seed = nmf_seed
  )
  saveRDS(geneNMF_programs, file.path(outdir, "03_geneNMF_programs.rds"))
  
  nmf_genes <- GeneNMF::getNMFgenes(
    nmf.res = geneNMF_programs,
    specificity.weight = 5,
    weight.explained = 0.5,
    max.genes = 100
  )
  saveRDS(nmf_genes, file.path(outdir, "04_nmf_genes.rds"))
  
  nmf_genes_long <- bind_rows(lapply(names(nmf_genes), function(prog) {
    data.frame(program = prog, gene = nmf_genes[[prog]], stringsAsFactors = FALSE)
  }))
  write.csv(nmf_genes_long, file.path(outdir, "04_nmf_genes_long.csv"), row.names = FALSE)
  
  #nMP_use <- max(4, min(8, floor(length(nmf_genes) / 4)))
  nMP_use <- 6
  geneNMF_metaprograms <- GeneNMF::getMetaPrograms(
    nmf.res = geneNMF_programs,
    nMP = nMP_use,
    specificity.weight = 5,
    weight.explained = 0.5,
    max.genes = 100,
    metric = "cosine",
    hclust.method = "ward.D2",
    min.confidence = 0.5,
    remove.empty = TRUE
  )
  saveRDS(geneNMF_metaprograms, file.path(outdir, "05_geneNMF_metaprograms.rds"))
  
  mp_genes <- geneNMF_metaprograms$metaprograms.genes
  mp_metrics <- geneNMF_metaprograms$metaprograms.metrics
  write.csv(mp_metrics, file.path(outdir, "05_metaprogram_metrics.csv"), row.names = FALSE)
  
  mp_genes_long <- bind_rows(lapply(names(mp_genes), function(mp) {
    data.frame(metaprogram = mp, gene = mp_genes[[mp]], stringsAsFactors = FALSE)
  }))
  write.csv(mp_genes_long, file.path(outdir, "05_metaprogram_genes_long.csv"), row.names = FALSE)
  
  DefaultAssay(seu_nm) <- assay_use
  pb_list <- AggregateExpression(object = seu_nm, assays = assay_use, group.by = sample_col, return.seurat = FALSE, verbose = FALSE)
  pb_counts <- as.matrix(pb_list[[assay_use]])
  
  colnames(pb_counts) <- standardize_sample_for_merge(colnames(pb_counts))
  keep_samples_std <- standardize_sample_for_merge(keep_samples)
  pb_counts <- pb_counts[, keep_samples_std[keep_samples_std %in% colnames(pb_counts)], drop = FALSE]
  
  rownames(pb_counts) <- toupper(rownames(pb_counts))
  pb_counts <- rowsum(pb_counts, group = rownames(pb_counts), reorder = FALSE)
  pb_counts <- as.matrix(pb_counts)
  pb_logcpm <- edgeR::cpm(pb_counts, log = TRUE, prior.count = 1)
  saveRDS(pb_counts, file.path(outdir, "07_pseudobulk_counts.rds"))
  saveRDS(pb_logcpm, file.path(outdir, "07_pseudobulk_logcpm.rds"))
  
  rank_data <- singscore::rankGenes(pb_logcpm, tiesMethod = "min")
  mp_singscore_long <- bind_rows(lapply(names(mp_genes), function(mp) {
    genes <- intersect(toupper(mp_genes[[mp]]), rownames(pb_logcpm))
    if (length(genes) < 3) return(NULL)
    ss <- singscore::simpleScore(rankData = rank_data, upSet = genes, knownDirection = TRUE)
    score_col <- if ("TotalScore" %in% colnames(ss)) "TotalScore" else if ("UpScore" %in% colnames(ss)) "UpScore" else colnames(ss)[1]
    tibble(sample = standardize_sample_for_merge(rownames(ss)), metaprogram = mp, singscore = ss[[score_col]], n_overlap = length(genes))
  }))
  write.csv(mp_singscore_long, file.path(outdir, "08_metaprogram_singscore_long.csv"), row.names = FALSE)
  
  mp_feature_wide <- mp_singscore_long %>%
    mutate(mp_key = sanitize_name(metaprogram)) %>%
    select(sample, mp_key, singscore, n_overlap) %>%
    pivot_wider(
      names_from = mp_key,
      values_from = c(singscore, n_overlap),
      names_glue = "{mp_key}__{.value}"
    )
  colnames(mp_feature_wide) <- paste0("denovo_gene_", colnames(mp_feature_wide))
  colnames(mp_feature_wide)[1] <- "sample"
  write.csv(mp_feature_wide, file.path(outdir, "08_metaprogram_feature_wide.csv"), row.names = FALSE)
  
  list(
    mp_singscore_long = mp_singscore_long,
    mp_feature_wide = mp_feature_wide,
    mp_genes = mp_genes,
    mp_metrics = mp_metrics
  )
}


merge_feature_blocks <- function(continuum_df = NULL, curated_feature_wide = NULL, nmf_feature_wide = NULL) {
  out <- NULL
  if (!is.null(continuum_df)) out <- continuum_df
  if (is.null(out) && !is.null(curated_feature_wide)) out <- curated_feature_wide
  if (is.null(out) && !is.null(nmf_feature_wide)) out <- nmf_feature_wide
  if (is.null(out)) stop("No feature blocks to merge.")
  
  if (!is.null(curated_feature_wide) && !identical(out, curated_feature_wide)) {
    out <- out %>% left_join(curated_feature_wide, by = "sample")
  }
  if (!is.null(nmf_feature_wide) && !identical(out, nmf_feature_wide)) {
    out <- out %>% left_join(nmf_feature_wide, by = "sample")
  }
  out
}

export_cellrank_input <- function(
    seu,
    outdir,
    assay_use = "RNA",
    cohort_col = "cohort_l2",
    target_cohort = "Lee_p1_base",
    sample_col = "sample",
    state_col = "manual.cluster_l2",
    keep_states = NULL,
    min_cells_per_patient = 20
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  md <- seu@meta.data
  
  keep_idx <- md[[cohort_col]] == target_cohort
  if (!is.null(keep_states)) keep_idx <- keep_idx & md[[state_col]] %in% keep_states
  keep_cells <- rownames(md)[keep_idx]
  seu_sub <- subset(seu, cells = keep_cells)
  
  patient_n <- as.data.table(seu_sub@meta.data)[, .N, by = sample_col]
  setnames(patient_n, "N", "n_cells_retained")
  eligible_patients <- patient_n[n_cells_retained >= min_cells_per_patient][[sample_col]]
  
  seu_sub <- subset(seu_sub, cells = rownames(seu_sub@meta.data)[seu_sub@meta.data[[sample_col]] %in% eligible_patients])
  patient_n2 <- as.data.table(seu_sub@meta.data)[, .N, by = sample_col]
  setnames(patient_n2, "N", "n_cells_retained")
  fwrite(patient_n2, file.path(outdir, "patient_counts_after_filter.csv"))
  
  counts <- GetAssayData(seu_sub, assay = assay_use, layer = "counts")
  gene_keep <- Matrix::rowSums(counts > 0) >= 3
  counts <- counts[gene_keep, ]
  
  obs <- seu_sub@meta.data
  obs$cell_id <- rownames(obs)
  obs <- obs[colnames(counts), , drop = FALSE]
  var <- data.table(gene_id = rownames(counts), gene_symbol = rownames(counts))
  
  writeMM(counts, file.path(outdir, "counts.mtx"))
  fwrite(as.data.table(obs), file.path(outdir, "obs.csv"))
  fwrite(var, file.path(outdir, "var.csv"))
  fwrite(data.table(cell_id = colnames(counts)), file.path(outdir, "barcodes.csv"))
  
  summary_dt <- data.table(
    n_cells = ncol(counts),
    n_genes = nrow(counts),
    n_patients = length(unique(obs[[sample_col]])),
    cohort = target_cohort
  )
  fwrite(summary_dt, file.path(outdir, "export_summary.csv"))
}

merge_cellrank_features <- function(feature_df, cellrank_patient_csv) {
  cellrank_df <- fread(cellrank_patient_csv)
  cellrank_df$sample <- standardize_sample_for_merge(cellrank_df$sample)
  patient_col <- "sample"
  
  keep_cols <- colnames(cellrank_df)[grepl("^cr_", colnames(cellrank_df)) | grepl("^qc_", colnames(cellrank_df))]
  cellrank_df <- cellrank_df %>%
    select(all_of(patient_col), all_of(keep_cols)) %>%
    distinct(.data[[patient_col]], .keep_all = TRUE)
  
  colnames(cellrank_df) <- paste0("cellrank_", colnames(cellrank_df))
  colnames(cellrank_df)[1] <- "sample"
  feature_df %>% left_join(cellrank_df, by = patient_col)
}

# =========================================================
# 3. LOAD CLEAN SUBSET
# =========================================================

seu_clean <- subset_clean
DefaultAssay(seu_clean) <- ASSAY_USE
check_clean_subset(seu_clean, SUBTYPE_COL, SAMPLE_COL, COHORT_COL, ASSAY_USE)

if (length(EXTRA_DROP_STATES) > 0) {
  seu_clean <- subset(seu_clean, subset = !(manual.cluster_l2 %in% EXTRA_DROP_STATES))
}

#saveRDS(seu_clean, file.path(OUTROOT, paste0("clean_input_confirmed_", sanitize_name(TARGET_ANALYSIS), ".rds")))

# =========================================================
# 4. PATIENT COMPOSITION
# =========================================================

comp_long <- make_patient_composition(
  seu = seu_clean,
  sample_col = SAMPLE_COL,
  cohort_col = COHORT_COL,
  response_col = RESPONSE_COL,
  celltype_col = SUBTYPE_COL
)

ggplot(comp_long, aes(x = manual.cluster_l2, y = n_cells)) +
  geom_boxplot(outlier.size = 0) +
  ggbeeswarm::geom_beeswarm(size = 1.2, alpha = 0.6, cex = 2.2) +
  facet_wrap(~manual.cluster_l2, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

comp_feature_wide <- comp_long %>%
  mutate(
    sample = standardize_sample_for_merge(sample),
    manual.cluster_l2 = as.character(manual.cluster_l2),
    prop = as.numeric(prop)
  ) %>%
  filter(!is.na(sample), sample != "", !is.na(manual.cluster_l2), manual.cluster_l2 != "") %>%
  group_by(sample, manual.cluster_l2) %>%
  summarise(prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = manual.cluster_l2,
    values_from = prop,
    values_fill = list(prop = 0),
    names_prefix = "composition_"
  )
# write.csv(comp_feature_wide, file.path(composition_outdir, "composition_feature_wide.csv"), row.names = FALSE)

# =========================================================
# 5. OPTIONAL CENTROID SCORE
# =========================================================
cat("\n[3/6] Centroid score...\n")
centroid_outdir <- file.path(OUTROOT, "centroid_from_clean_subset")
dir.create(centroid_outdir, recursive = TRUE, showWarnings = FALSE)

ElbowPlot(subset_clean)#20260327_monocyte_clean_pca_elbow
DimHeatmap(subset_clean, dims = 1:12, cells = 500, balanced = TRUE) #8x8 #20260327_monocyte_clean_pca

continuum_df <- NULL
if (!is.na(ANCHOR_STATE_LOW) && !is.na(ANCHOR_STATE_HIGH)) {
  centroid_res <- calc_two_anchor_score(
    seu = subset(seu_clean, subset = cohort_l2 == TARGET_COHORT),
    sample_col = SAMPLE_COL,
    celltype_col = SUBTYPE_COL,
    low_label = ANCHOR_STATE_LOW,
    high_label = ANCHOR_STATE_HIGH,
    assay_use = ASSAY_USE,
    rerun_pca = TRUE,
    dims_use = c(1,4,7,8,9),
    npcs = 20,
    min_cells_per_patient = MIN_CELLS_PROP,
    prefix = CENTROID_PREFIX
  )
  saveRDS(centroid_res, file.path(centroid_outdir, paste0(sanitize_name(CENTROID_PREFIX), "_res.rds")))
  
  continuum_df <- centroid_res$patient_summary %>%
    rename(sample = .sample) %>%
    mutate(sample = standardize_sample_for_merge(sample))
  colnames(continuum_df) <- paste0("Centroid_", colnames(continuum_df))
  colnames(continuum_df)[1] <- "sample"
  write.csv(continuum_df, file.path(centroid_outdir, "centroid_feature_wide.csv"), row.names = FALSE)
} else {
  cat("Skipping centroid score because ANCHOR_STATE_LOW/HIGH are NA.\n")
}

# =========================================================
# 6. CURATED PSEUDOBULK + NMF
# =========================================================
cat("\n[4/6] Curated pseudobulk scores...\n")
curated_outdir <- file.path(OUTROOT, "curated_program_scores_from_clean_subset")
curated_res <- run_curated_pseudobulk_scores(
  seu_use = seu_clean,
  outdir = curated_outdir,
  sample_col = SAMPLE_COL,
  cohort_col = COHORT_COL,
  target_cohort = TARGET_COHORT,
  assay_use = ASSAY_USE,
  species_use = "Homo sapiens",
  gene_sets = CUSTOM_GENE_SETS,
  hallmark_targets = MSIGDB_TARGETS,
  min_cells_per_sample = MIN_CELLS_PB,
  min_nonzero_cells_per_sample = MIN_NONZERO_PER_SMP,
  min_samples_passing_detection = MIN_SAMPLES_DETECT,
  min_geneset_overlap = MIN_GENESET_OVERLAP,
  prefix_name = sanitize_name(TARGET_ANALYSIS)
)

curated_feature_wide <- curated_res$feature_long %>%
  dplyr::mutate(
    sample = as.character(sample),
    gs_key = sanitize_name(gs_name)
  ) %>%
  dplyr::filter(
    !is.na(sample), sample != "",
    !is.na(gs_key), gs_key != ""
  ) %>%
  dplyr::group_by(sample, gs_key) %>%
  dplyr::summarise(
    ssgsea2_es = if (all(is.na(ssgsea2_es))) NA_real_ else mean(ssgsea2_es, na.rm = TRUE),
    singscore  = if (all(is.na(singscore)))  NA_real_ else mean(singscore,  na.rm = TRUE),
    pc1        = if (all(is.na(pc1)))        NA_real_ else mean(pc1,        na.rm = TRUE),
    pc2        = if (all(is.na(pc2)))        NA_real_ else mean(pc2,        na.rm = TRUE),
    eigengene  = if (all(is.na(eigengene)))  NA_real_ else mean(eigengene,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = gs_key,
    values_from = c(ssgsea2_es, singscore, pc1, pc2, eigengene),
    names_glue = "{gs_key}__{.value}"
  ) %>%
  tibble::as_tibble()

colnames(curated_feature_wide) <- paste0("curated_gene_", colnames(curated_feature_wide))
colnames(curated_feature_wide)[1] <- "sample"

curated_res$curated_feature_wide <- curated_feature_wide

dim(curated_res$curated_feature_wide)
head(curated_res$curated_feature_wide[, 1:min(8, ncol(curated_res$curated_feature_wide))])
curated_res$curated_feature_wide$sample <-gsub('_','-',curated_res$curated_feature_wide$sample)
cat("\n[5/6] GeneNMF...\n")
nmf_outdir <- file.path(OUTROOT, "geneNMF_from_clean_subset")
nmf_res <- run_gene_nmf_features(
  seu_use = seu_clean,
  outdir = nmf_outdir,
  sample_col = SAMPLE_COL,
  cohort_col = COHORT_COL,
  target_cohort = TARGET_COHORT,
  assay_use = ASSAY_USE,
  min_cells_per_sample = MIN_CELLS_PB,
  k_vec = 4:6,
  nfeatures_hvg = 2000,
  nmf_seed = 123,
  gene_blocklist_patterns = GENE_BLOCKLIST_PATTERNS
)


cat("\n[6/6] Merging features...\n")
final_feature_df <- merge_feature_blocks(
  continuum_df = if (!is.null(continuum_df)) continuum_df else comp_feature_wide,
  curated_feature_wide = curated_res$curated_feature_wide,
  nmf_feature_wide = nmf_res$mp_feature_wide
)
#########
# 2. composition
########
# naive_cytotoxic <- subset(CD8T_subset, subset = manual.cluster_l2 %in% c("cytotoxic_CD8","naive_CD8"))
# naive_cytotoxic <- subset(naive_cytotoxic, subset = cohort_l2 == 'Lee_p1_base')

make_patient_ilr <- function(
    seu,
    sample_col,
    celltype_col,
    pseudocount = 0.5,
    celltype_levels = NULL,
    min_cells_per_type = 10
) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))
  
  ilr_from_prop <- function(prop_mat) {
    x <- as.matrix(prop_mat)
    if (ncol(x) < 2) stop("ILR requires at least 2 cell types.")
    
    x <- x / rowSums(x)
    logx <- log(x)
    clr <- logx - rowMeans(logx)
    
    V <- as.matrix(stats::contr.helmert(ncol(x)))
    V <- apply(V, 2, function(v) v / sqrt(sum(v^2)))
    
    ilr <- clr %*% V
    ilr <- as.data.frame(ilr)
    colnames(ilr) <- paste0("ilr_", seq_len(ncol(ilr)))
    ilr
  }
  
  meta_df <- seu@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")
  
  # sample_col이 비어 있으면 donor_id로 보정
  if ("donor_id" %in% colnames(meta_df)) {
    meta_df[[sample_col]] <- ifelse(
      is.na(meta_df[[sample_col]]) | meta_df[[sample_col]] == "",
      meta_df$donor_id,
      meta_df[[sample_col]]
    )
  }
  
  # sample x celltype count long table
  comp_long_full <- meta_df %>%
    dplyr::filter(!is.na(.data[[sample_col]]), .data[[sample_col]] != "") %>%
    dplyr::count(
      .data[[sample_col]],
      .data[[celltype_col]],
      name = "n_cells"
    ) %>%
    dplyr::rename(
      sample = dplyr::all_of(sample_col),
      manual.cluster_l2 = dplyr::all_of(celltype_col)
    )
  
  if (nrow(comp_long_full) == 0) {
    stop("No valid cells found after filtering missing sample IDs.")
  }
  
  # celltype 순서 고정
  if (is.null(celltype_levels)) {
    celltype_levels <- sort(unique(comp_long_full$manual.cluster_l2))
  }
  
  # sample x celltype wide count matrix
  # sample별로 실제 존재하는 조합 안에서만 subtype을 채움
  comp_wide_counts_full <- comp_long_full %>%
    dplyr::mutate(
      manual.cluster_l2 = factor(manual.cluster_l2, levels = celltype_levels)
    ) %>%
    tidyr::complete(
      tidyr::nesting(sample),
      manual.cluster_l2,
      fill = list(n_cells = 0)
    ) %>%
    tidyr::pivot_wider(
      names_from = manual.cluster_l2,
      values_from = n_cells,
      values_fill = 0
    ) %>%
    dplyr::arrange(sample)
  
  # 모든 지정 subtype이 min_cells_per_type 이상인 sample만 유지
  keep_sample_df <- comp_wide_counts_full %>%
    dplyr::mutate(
      keep_for_ilr = dplyr::if_all(
        dplyr::all_of(celltype_levels),
        ~ .x >= min_cells_per_type
      )
    ) %>%
    dplyr::select(sample, keep_for_ilr)
  
  comp_wide_counts <- comp_wide_counts_full %>%
    dplyr::semi_join(
      keep_sample_df %>% dplyr::filter(keep_for_ilr),
      by = "sample"
    )
  
  kept_samples <- comp_wide_counts %>% dplyr::pull(sample)
  excluded_samples <- setdiff(comp_wide_counts_full$sample, kept_samples)
  
  if (nrow(comp_wide_counts) == 0) {
    stop("No samples passed the min_cells_per_type filter. Lower the threshold or reduce the number of cell types.")
  }
  
  # filtered long table + proportion
  comp_long <- comp_long_full %>%
    dplyr::filter(sample %in% kept_samples) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      total_cells = sum(n_cells),
      prop = n_cells / total_cells
    ) %>%
    dplyr::ungroup()
  
  # count matrix
  count_mat <- comp_wide_counts %>%
    dplyr::select(dplyr::all_of(celltype_levels)) %>%
    as.matrix()
  rownames(count_mat) <- comp_wide_counts$sample
  
  # pseudocount 추가 후 proportion 계산
  count_mat_pc <- count_mat + pseudocount
  prop_mat_pc <- count_mat_pc / rowSums(count_mat_pc)
  
  # ILR
  ilr_df <- ilr_from_prop(prop_mat_pc)
  ilr_df <- dplyr::bind_cols(
    comp_wide_counts %>% dplyr::select(sample),
    ilr_df
  )
  
  # pseudocount 반영 proportion wide
  prop_wide_pseudocount <- as.data.frame(prop_mat_pc)
  prop_wide_pseudocount <- dplyr::bind_cols(
    comp_wide_counts %>% dplyr::select(sample),
    prop_wide_pseudocount
  )
  
  list(
    comp_long = comp_long,
    comp_long_full = comp_long_full,
    comp_wide_counts = comp_wide_counts,
    comp_wide_counts_full = comp_wide_counts_full,
    prop_wide_pseudocount = prop_wide_pseudocount,
    ilr_df = ilr_df,
    keep_sample_df = keep_sample_df,
    kept_samples = kept_samples,
    excluded_samples = excluded_samples,
    celltype_levels = celltype_levels,
    pseudocount = pseudocount,
    min_cells_per_type = min_cells_per_type
  )
}
SUBTYPE_COL = "manual.cluster_l2"

res_ilr <- make_patient_ilr(
  seu = subset_clean,
  sample_col = SAMPLE_COL,
  celltype_col = SUBTYPE_COL,
  pseudocount = 0.5,
  celltype_levels = c(
    "activated monocyte","CD16 monotype"
  ),
  min_cells_per_type = 10
)

ilr_df <- res_ilr$ilr_df
comp_long <- res_ilr$comp_long
res_ilr$keep_sample_df
ilr_df$sample <- gsub('_','-',ilr_df$sample)
#########
# final_feature_df <- readRDS('data/20260309_pilot/results/version2/NK/final_feature_df_nk_from_clean_subset.rds')
# final_feature_df <- as.data.frame(final_feature_df)
# final_feature_df <- final_feature_df[,-c(101,102)]
# final_feature_df <- merge(final_feature_df, ilr_df,by='sample',all.x=T)
#
# always add composition block
if (!all(colnames(comp_feature_wide) %in% colnames(final_feature_df))) {
  final_feature_df <- final_feature_df %>% left_join(comp_feature_wide, by = "sample")
}
final_feature_df <- final_feature_df %>% distinct(sample, .keep_all = TRUE)
write.csv(final_feature_df, file.path(OUTROOT, paste0("final_feature_df_", sanitize_name(TARGET_ANALYSIS), "_from_clean_subset.csv")), row.names = FALSE)
saveRDS(final_feature_df, file.path(OUTROOT, paste0("final_feature_df_", sanitize_name(TARGET_ANALYSIS), "_from_clean_subset.rds")))

if (MERGE_WITH_EXISTING_FINAL_FEATURE && nzchar(EXISTING_FINAL_FEATURE_CSV) && file.exists(EXISTING_FINAL_FEATURE_CSV)) {
  existing_df <- fread(EXISTING_FINAL_FEATURE_CSV) %>% as.data.frame()
  existing_df$sample <- standardize_sample_for_merge(existing_df$sample)
  merged_df <- existing_df %>% left_join(final_feature_df, by = "sample")
  write.csv(merged_df, file.path(OUTROOT, paste0("final_feature_df_merged_into_existing_", sanitize_name(TARGET_ANALYSIS), ".csv")), row.names = FALSE)
  saveRDS(merged_df, file.path(OUTROOT, paste0("final_feature_df_merged_into_existing_", sanitize_name(TARGET_ANALYSIS), ".rds")))
}

if (DO_CELLRANK_EXPORT) {
  cellrank_outdir <- file.path(OUTROOT, "CellRank_input_from_clean_subset")
  export_cellrank_input(
    seu = seu_clean,
    outdir = cellrank_outdir,
    assay_use = ASSAY_USE,
    cohort_col = COHORT_COL,
    target_cohort = TARGET_COHORT,
    sample_col = SAMPLE_COL,
    state_col = SUBTYPE_COL,
    keep_states = KEEP_STATES_FOR_CELLRANK,
    min_cells_per_patient = MIN_CELLS_PB
  )
}

if (nzchar(CELLRANK_PATIENT_CSV) && file.exists(CELLRANK_PATIENT_CSV)) {
  final_feature_df_cellrank <- merge_cellrank_features(final_feature_df, CELLRANK_PATIENT_CSV)
  write.csv(final_feature_df_cellrank, file.path(OUTROOT, paste0("final_feature_df_with_cellrank_", sanitize_name(TARGET_ANALYSIS), ".csv")), row.names = FALSE)
  saveRDS(final_feature_df_cellrank, file.path(OUTROOT, paste0("final_feature_df_with_cellrank_", sanitize_name(TARGET_ANALYSIS), ".rds")))
}

final_feature_df_cellrank_filt <- final_feature_df_cellrank[,!(colnames(final_feature_df_cellrank) %in% c('denovo_gene_mp2__n_overlap', 'denovo_gene_mp3__n_overlap', 'denovo_gene_mp4__n_overlap', 'denovo_gene_mp5__n_overlap', 'denovo_gene_mp6__n_overlap',  'denovo_gene_mp7__n_overlap','Centroid_n_score_cells',  'Centroid_keep_for_analysis', 'cellrank_cr_b3_priming_mean', 'cellrank_cr_b3_priming_q75', 'cellrank_cr_b3_n_cells', 'cellrank_qc_b3_subset_frac_CD16 monotype','cellrank_qc_b3_subset_frac_activated monocyte'))]
OUT_CSV <- "data/20260309_pilot/results/version2/Monocyte/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv"
fwrite(as.data.table(final_feature_df_cellrank_filt), OUT_CSV)

