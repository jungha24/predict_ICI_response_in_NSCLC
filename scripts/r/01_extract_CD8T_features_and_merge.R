# 2026.03.30
# make feature library in CD4T
# noncanonical subset 따로 떼서 다음 분석에 이용하기
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
DimPlot(seu_integrated, reduction = "umap", group.by='manual.cluster')
table(seu_integrated@meta.data[seu_integrated@meta.data$manual.cluster=='naive/TCM CD4+T',]$predicted.celltype.l1)
table(seu_integrated@meta.data[seu_integrated@meta.data$manual.cluster %in% c('navie CD8+T','effector CD8+T'),]$predicted.celltype.l1)
DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l1')

## subclustering
CD8T_subset<- subset(seu_integrated, idents = c('effector CD8+T','navie CD8+T'))
DefaultAssay(CD8T_subset) <- 'RNA'
CD8T_subset$source_obj <- ifelse(grepl("^DIS_", colnames(CD8T_subset)), "disease", "normal")
CD8T_list <- SplitObject(CD8T_subset, split.by = "source_obj")
names(CD8T_list)

CD8T_list <- lapply(CD8T_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) #<<<
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000) #<<<
  return(x)
})

features <- SelectIntegrationFeatures(object.list = CD8T_list, nfeatures = 3000) #<<<

CD8T_list <- lapply(CD8T_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

anchors <- FindIntegrationAnchors(
  object.list = CD8T_list,
  anchor.features = features,
  reduction = "rpca", #<<<
  dims = 1:30 #<<<
)

CD8T_subset <- IntegrateData(
  anchorset = anchors,
  dims = 1:30 #<<<
)

CD8T_subset <- ScaleData(CD8T_subset, verbose = FALSE)
CD8T_subset <- RunPCA(CD8T_subset, npcs = 50, verbose = TRUE)
ElbowPlot(CD8T_subset) #3x4 #20260324_B_subcluster_PCA_elbowplot

CD8T_subset <- FindNeighbors(CD8T_subset, dims = 1:8) #<<< 
CD8T_subset <- FindClusters(CD8T_subset, resolution = 0.1) #<<<
CD8T_subset <- RunUMAP(CD8T_subset, dims = 1:8) #<<<
DimPlot(CD8T_subset, reduction='umap',label=T)#3x4
DimPlot(CD8T_subset, reduction='umap', group.by = 'predicted.celltype.l2',label=T)#3x4
markers <- FindAllMarkers(CD8T_subset, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

DimPlot(CD8T_subset, group.by='cohort_l2', split.by = 'cohort_l2')

CD8T_subset$manual.cluster_l2 <- ifelse(CD8T_subset$seurat_clusters==0,'cytotoxic_CD8',
                                        ifelse(CD8T_subset$seurat_clusters==1, 'naive_CD8',
                                               ifelse(CD8T_subset$seurat_clusters==2,'unconventional_T',
                                                      ifelse(CD8T_subset$seurat_clusters==3,'MAIT','cycling_T'
                                                             ))))

cytotoxic_marker <- c('FGFBP2','GZMH','LGALS1','ITGB1')
naive_marker <- c('CCR7','LEF1','TCF7','SELL')

DimPlot(CD8T_subset, label=T, group.by='manual.cluster_l2')#20260331_CD8T_subcluster_umap
saveRDS(CD8T_subset, file='data/20260309_pilot/results/version2/CD8T/CD8T_subset.rds')
#####
# =========================================================
# 0. USER CONFIG
# =========================================================
TARGET_ANALYSIS <- "CD8T"   # "Monocyte" or "NK"

# optional: reload clinical metadata from original meta file
ADD_CLINICAL_META <- TRUE
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

# optional merge with existing patient feature matrix
MERGE_WITH_EXISTING_FINAL_FEATURE <- FALSE
EXISTING_FINAL_FEATURE_CSV <- ""

# optional CellRank export / merge
DO_CELLRANK_EXPORT <- TRUE
KEEP_STATES_FOR_CELLRANK <- NULL  # NULL = use all states present in the clean subset
CELLRANK_PATIENT_CSV <- ""
# if provided and exists, merge CellRank patient features

# output root
OUTROOT <- file.path("data/20260309_pilot/results/version2", TARGET_ANALYSIS)
dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)
#########
# 1. continuum feature
########
naive_cytotoxic <- subset(CD8T_subset, subset = manual.cluster_l2 %in% c("cytotoxic_CD8","naive_CD8"))
naive_cytotoxic <- subset(naive_cytotoxic, subset = cohort_l2 == 'Lee_p1_base')

DefaultAssay(naive_cytotoxic) <- 'RNA'
naive_cytotoxic<- JoinLayers(naive_cytotoxic)

naive_cytotoxic <- NormalizeData(naive_cytotoxic, normalization.method = "LogNormalize", scale.factor = 10000) 
naive_cytotoxic <- FindVariableFeatures(naive_cytotoxic, selection.method = "vst", nfeatures = 3000) #<<<

naive_cytotoxic <- ScaleData(naive_cytotoxic, verbose = FALSE)
naive_cytotoxic <- RunPCA(naive_cytotoxic, npcs = 50, verbose = TRUE)

DimPlot(naive_cytotoxic, reduction='pca', group.by='manual.cluster_l2')
ElbowPlot(naive_cytotoxic) #3x4 #20260330_CD8T_naive_cytotoxic_PCA_elbow
DimHeatmap(naive_cytotoxic, dims = 1:10, cells = 500, balanced = TRUE) #8x8 #20260330_CD8T_naive_cytotoxic_PCA_dimheat
# PC: 1, 4
ANCHOR_STATE_LOW  <- 'naive_CD8'
ANCHOR_STATE_HIGH <- 'cytotoxic_CD8'
CENTROID_PREFIX   <-"cd8t_naive_to_cytotoxic"

centroid_outdir <- file.path(OUTROOT, "centroid_from_clean_subset")
dir.create(centroid_outdir, recursive = TRUE, showWarnings = FALSE)

continuum_df <- NULL
if (!is.na(ANCHOR_STATE_LOW) && !is.na(ANCHOR_STATE_HIGH)) {
  continuum_res <- calc_two_anchor_score(
    seu = naive_cytotoxic,
    sample_col = "sample",
    celltype_col = "manual.cluster_l2",
    cohort_col = "cohort_l2",
    target_cohort = "Lee_p1_base",
    low_label = ANCHOR_STATE_LOW,
    high_label = ANCHOR_STATE_HIGH,
    assay_use = "RNA",
    rerun_pca = TRUE,
    dims_use = c(1,4),
    npcs = 20,
    min_cells_per_sample_input = 10,
    min_cells_per_patient = 10,
    prefix = CENTROID_PREFIX
  )
  saveRDS(continuum_res, file.path(centroid_outdir, paste0(sanitize_name(CENTROID_PREFIX), "_res.rds")))
  
  continuum_df <- continuum_res$patient_summary %>%
    rename(sample = .sample) %>%
    mutate(sample = standardize_sample_for_merge(sample))
  colnames(continuum_df) <- paste0("Centroid_", colnames(continuum_df))
  colnames(continuum_df)[1] <- "sample"
  write.csv(continuum_df, file.path(centroid_outdir, "centroid_feature_wide.csv"), row.names = FALSE)
} else {
  cat("Skipping centroid score because ANCHOR_STATE_LOW/HIGH are NA.\n")
}


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
res_ilr <- make_patient_ilr(
  seu = naive_cytotoxic,
  sample_col = SAMPLE_COL,
  celltype_col = SUBTYPE_COL,
  pseudocount = 0.5,
  celltype_levels = c(
    "cytotoxic_CD8","naive_CD8"
  ),
  min_cells_per_type = 10
)

ilr_df <- res_ilr$ilr_df
comp_long <- res_ilr$comp_long
res_ilr$keep_sample_df
ilr_df$sample <- gsub('_','-',ilr_df$sample)

#########
# 3. cellrank
########


if (DO_CELLRANK_EXPORT) {
  cellrank_outdir <- file.path(OUTROOT, "CellRank_input_from_clean_subset")
  export_cellrank_input(
    seu = naive_cytotoxic,
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

cytotoxic_marker <- c('FGFBP2','GZMH','LGALS1','ITGB1','NKG7','GNLY','PRF1','CCL5')
naive_marker <- c('CCR7','LEF1','TCF7','SELL','IL7R')

VlnPlot(naive_cytotoxic, features=cytotoxic_marker,pt.size=0)
VlnPlot(naive_cytotoxic, features=naive_marker,pt.size=0) # cluster 1= resting


#########
# 4. curated gene set
########

get_analysis_defaults <- function(target_analysis) {
  norm_symbol <- function(x) {
    unique(toupper(trimws(as.character(x))))
  }
  
  if (target_analysis == "CD8T") {
    
    list(
      gene_blocklist_patterns = c(
        "^MT-", "^RPS", "^RPL", "^HSP",
        "^IG[HKL]",
        "^TRA", "^TRB", "^TRD", "^TRG"
      ),
      custom_gene_sets = list(
        CD8_NAIVE_CENTRAL_MEMORY = c(
          "CCR7", "TCF7", "LEF1", "SELL", "IL7R",
          "LTB", "MAL", "SATB1", "BACH2", "TXK",
          "NOSIP", "LTB", "TMIGD2"
        ),
        
        CD8_GZMK_EFFECTOR_MEMORY = c(
          "GZMK", "CCL5", "CXCR3", "IL7R", "AQP3",
          "KLRB1", "LTB", "HOPX", "CXCR6", "DUSP2",
          "EOMES", "IFITM3"
        ),
        
        CD8_CYTOTOXIC_EFFECTOR = c(
          "NKG7", "PRF1", "GZMB", "GZMH", "GNLY",
          "CCL5", "CTSW", "FGFBP2", "CX3CR1", "KLRD1",
          "FCGR3A", "HOPX", "TBX21"
        ),
        
        CD8_TERMINAL_TEMRA_NKLIKE = c(
          "FGFBP2", "GZMH", "GNLY", "PRF1", "CCL5",
          "CX3CR1", "FCGR3A", "KLRD1", "ADGRG1", "CTSW",
          "ZEB2", "TBX21", "IFITM3"
        ),
        
        CD8_TCR_ACTIVATION_IMMEDIATE = c(
          "CD69", "TNFRSF9", "TNFRSF4",
          "NR4A1", "NR4A2", "NR4A3",
          "FOS", "FOSB", "JUN", "JUNB",
          "EGR1", "EGR2", "NFKBIA", "DUSP1",
          "ZFP36", "PPP1R15A"
        ),
        
        CD8_IFN_RESPONSE = c(
          "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6",
          "IFI44", "IFI44L", "MX1", "MX2", "OAS1",
          "OASL", "STAT1", "IRF7", "GBP1", "CXCL10"
        ),
        
        CD8_DYSFUNCTION_EXHAUSTION = c(
          "PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4",
          "TOX", "TOX2", "ENTPD1", "CXCL13", "BATF",
          "PRDM1", "MAF", "LAYN"
        ),
        
        CD8_PROLIFERATION_CYCLING = c(
          "MKI67", "TYMS", "PCLAF", "RRM2", "TK1",
          "STMN1", "HMGB2", "BIRC5", "MCM5", "MCM6",
          "PTTG1", "TOP2A"
        ),
        
        CD8_TISSUE_RESIDENT_LIKE = c(
          "ITGAE", "CXCR6", "ZNF683", "XCL1", "XCL2",
          "HOBIT", "RUNX3", "PDCD1", "CD69", "CRTAM"
        )
      ),
      msigdb_targets = c(
        "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_MTORC1_SIGNALING",
        "HALLMARK_GLYCOLYSIS",
        "HALLMARK_HYPOXIA",
        "HALLMARK_APOPTOSIS",
        "GOBP_T_CELL_ACTIVATION",
        "GOBP_ALPHA_BETA_T_CELL_ACTIVATION",
        "GOBP_REGULATION_OF_T_CELL_ACTIVATION",
        "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
        "GOBP_LYMPHOCYTE_MIGRATION",
        "GOBP_POSITIVE_REGULATION_OF_CELL_KILLING",
        "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
        "GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY",
        "REACTOME_TCR_SIGNALING",
        "REACTOME_DOWNSTREAM_TCR_SIGNALING",
        "REACTOME_CO_STIMULATION_BY_CD28",
        "REACTOME_INTERLEUKIN_7_SIGNALING",
        "REACTOME_SIGNALING_BY_INTERLEUKINS",
        "REACTOME_PD_1_SIGNALING"
      )
    )
    
  } else {
    stop("TARGET_ANALYSIS must be 'CD4T' or 'CD8T'.")
  }
}
cfg <- get_analysis_defaults(TARGET_ANALYSIS)
GENE_BLOCKLIST_PATTERNS <- cfg$gene_blocklist_patterns
CUSTOM_GENE_SETS <- cfg$custom_gene_sets
MSIGDB_TARGETS <- cfg$msigdb_targets

curated_outdir <- file.path(OUTROOT, "curated_program_scores_from_clean_subset")
curated_res <- run_curated_pseudobulk_scores(
  seu_use = naive_cytotoxic,
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
curated_feature_wide$sample <- gsub('_','-',curated_feature_wide$sample)
curated_res$curated_feature_wide <- curated_feature_wide

######
# NMF
#####
k_vec <- 4:7
nfeatures_hvg <- 2000
nmf_seed <- 123
sample_col <- "sample"
cohort_col <- "cohort_l2"
target_cohort <- "Lee_p1_base"
assay_use <- "RNA"
gene_blocklist_patterns <- c("^MT-","^RPS", "^RPL","^HSP","^IG[HKL]","^TRA", "^TRB", "^TRD", "^TRG")
min_cells_per_sample <- 20
meta0 <- naive_cytotoxic@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

## =========================================================
## 3. sample filtering
## =========================================================
meta1 <- naive_cytotoxic@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

keep_cells_non_na <- meta1$cell[!is.na(meta1[[sample_col]]) & meta1[[sample_col]] != ""]
seu_non_na <- subset(naive_cytotoxic, cells = keep_cells_non_na)

meta2 <- seu_non_na@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

sample_meta2 <- meta2 %>%
  transmute(sample = as.character(.data[[sample_col]])) %>%
  count(sample, name = "n_cells") %>%
  arrange(desc(n_cells))

keep_samples <- sample_meta2 %>%
  filter(n_cells >= min_cells_per_sample) %>%
  pull(sample)

if (length(keep_samples) < 2) {
  stop("Need at least 2 samples for robust multi-sample GeneNMF consensus analysis.")
}

keep_cells <- meta2$cell[meta2[[sample_col]] %in% keep_samples]
seu_nm <- subset(seu_non_na, cells = keep_cells)

obj_list <- SplitObject(seu_nm, split.by = sample_col)

sample_sizes <- sapply(obj_list, ncol)
print(sample_sizes)

obj_list <- obj_list[sample_sizes >= min_cells_per_sample]
print(names(obj_list))

all_genes <- rownames(obj_list[[1]])
genes_blocklist <- unique(unlist(lapply(gene_blocklist_patterns, function(p) {
  grep(p, all_genes, value = TRUE, ignore.case = FALSE)
})))

print(length(genes_blocklist))

obj_list <- lapply(obj_list, function(x) {
  DefaultAssay(x) <- assay_use
  print(Layers(x[[assay_use]]))
  x <- JoinLayers(x, assay = assay_use)
  print(Layers(x[[assay_use]]))
  if (!"data" %in% Layers(x[[assay_use]])) {
    x <- NormalizeData(x, assay = assay_use, verbose = FALSE)
  }
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
hvg_summary <- data.frame(
  sample = names(hvg_list),
  n_hvg = sapply(hvg_list, length)
)

hvg_union <- sort(unique(unlist(hvg_list)))

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

nmf_genes <- GeneNMF::getNMFgenes(
  nmf.res = geneNMF_programs,
  specificity.weight = 5,
  weight.explained = 0.5,
  max.genes = 100
)

geneNMF_metaprograms <- GeneNMF::getMetaPrograms(
  nmf.res = geneNMF_programs,
  nMP = 5,                  # 일단 시작값; 필요시 4~10 비교
  specificity.weight = 5,
  weight.explained = 0.5,
  max.genes = 100,
  metric = "cosine",
  hclust.method = "ward.D2",
  min.confidence = 0.5,
  remove.empty = TRUE
)

mp_genes <- geneNMF_metaprograms$metaprograms.genes
mp_metrics <- geneNMF_metaprograms$metaprograms.metrics

anno_colors <- brewer.pal(n=5, name="Paired")
names(anno_colors) <- names(geneNMF_metaprograms$metaprograms.genes)
ph <- plotMetaPrograms(geneNMF_metaprograms, annotation_colors = anno_colors)
ph #8x9 #20260331_CD8T_NMF_metaprogam_heatmap

nmf_outdir <- file.path(OUTROOT, "geneNMF_from_clean_subset")
nmf_res <- run_gene_nmf_features(
  seu_use = naive_cytotoxic,
  outdir = nmf_outdir,
  sample_col = SAMPLE_COL,
  cohort_col = COHORT_COL,
  target_cohort = TARGET_COHORT,
  assay_use = ASSAY_USE,
  min_cells_per_sample = MIN_CELLS_PB,
  k_vec = 4:7,
  nfeatures_hvg = 2000,
  nmf_seed = 123,
  gene_blocklist_patterns = GENE_BLOCKLIST_PATTERNS
)
geneNMF_metaprograms <- readRDS('data/20260309_pilot/results/version2/CD4T/geneNMF_from_clean_subset/05_geneNMF_metaprograms.rds')
mp_genes <- geneNMF_metaprograms$metaprograms.genes
mp_metrics <- geneNMF_metaprograms$metaprograms.metrics

anno_colors <- brewer.pal(n=5, name="Paired")
names(anno_colors) <- names(geneNMF_metaprograms$metaprograms.genes)
ph <- plotMetaPrograms(geneNMF_metaprograms, annotation_colors = anno_colors)
ph 
##########
# merge
#######
final_feature_df <- merge_feature_blocks(
  continuum_df = if (!is.null(continuum_df)) continuum_df else comp_feature_wide,
  curated_feature_wide = curated_res$curated_feature_wide,
  nmf_feature_wide = nmf_res$mp_feature_wide
)

if (!all(colnames(ilr_df) %in% colnames(final_feature_df))) {
  final_feature_df <- final_feature_df %>% left_join(ilr_df, by = "sample")
}
final_feature_df <- final_feature_df %>% distinct(sample, .keep_all = TRUE)
write.csv(final_feature_df, file.path(OUTROOT, paste0("final_feature_df_", sanitize_name(TARGET_ANALYSIS), "_from_clean_subset.csv")), row.names = FALSE)
saveRDS(final_feature_df, file.path(OUTROOT, paste0("final_feature_df_", sanitize_name(TARGET_ANALYSIS), "_from_clean_subset.rds")))
final_feature_df <-readRDS('data/20260309_pilot/results/version2/CD8T/final_feature_df_cd8t_from_clean_subset.rds')
CELLRANK_PATIENT_CSV = 'data/20260309_pilot/results/version2/CD8T/CellRank_output_from_clean_subset/cellrank_b3_patient_features.csv'
if (nzchar(CELLRANK_PATIENT_CSV) && file.exists(CELLRANK_PATIENT_CSV)) {
  final_feature_df_cellrank <- merge_cellrank_features(final_feature_df, CELLRANK_PATIENT_CSV)
  write.csv(final_feature_df_cellrank, file.path(OUTROOT, paste0("final_feature_df_with_cellrank_", sanitize_name(TARGET_ANALYSIS), ".csv")), row.names = FALSE)
  saveRDS(final_feature_df_cellrank, file.path(OUTROOT, paste0("final_feature_df_with_cellrank_", sanitize_name(TARGET_ANALYSIS), ".rds")))
}

final_feature_df_cellrank_filt <- final_feature_df_cellrank[,!(colnames(final_feature_df_cellrank) %in% c('denovo_gene_mp1__singscore','denovo_gene_mp1__n_overlap','denovo_gene_mp2__n_overlap', 'denovo_gene_mp3__n_overlap', 'denovo_gene_mp4__n_overlap', 'denovo_gene_mp5__n_overlap', 'Centroid_n_score_cells',  'Centroid_keep_for_analysis', 'cellrank_cr_b3_priming_mean', 'cellrank_cr_b3_priming_q75','Centroid_cd4t_ier_to_active_n_score_cells','cellrank_qc_b3_subset_frac_cytotoxic_CD8','cellrank_qc_b3_subset_frac_naive_CD8'))]
OUT_CSV <- "data/20260309_pilot/results/version2/CD8T/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv"
fwrite(as.data.table(final_feature_df_cellrank_filt), OUT_CSV)
