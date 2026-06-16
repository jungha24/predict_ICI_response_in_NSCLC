#

load('data/20260309_pilot/results/20260310_integrated.RData')
selected_patients <- fread('data/20260309_pilot/results/20260313_patients_selected.txt',header=F)
meta <- fread('data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt')
## meta
idx <- match(seu_integrated@meta.data$sample, meta$sample_id)
seu_integrated$binarized_response <- meta$`Binarized response`[idx]
seu_integrated$`PD.L1_TPS` <- meta$`PD-L1_TPS_mod`[idx]

cd4t <- readRDS('data/20260309_pilot/results/version2/CD4T/CD4T_subset.rds')
DimPlot(cd4t, reduction='umap',group.by='manual.cluster_l2',label=T)
cd8t <- readRDS('data/20260309_pilot/results/version2/CD8T/CD8T_subset.rds')
DimPlot(cd8t, reduction='umap',group.by='manual.cluster_l2',label=T)

keep_cells_cd4t <- colnames(subset(cd4t, subset = manual.cluster_l2 == 'unconventional_T'))
keep_cells_cd8t <- colnames(subset(cd8t, subset = manual.cluster_l2 %in% c('unconventional_T','MAIT')))
nonconventional_T <- subset(seu_integrated, cells = c(keep_cells_cd4t,keep_cells_cd8t))

DefaultAssay(nonconventional_T) <- 'RNA'
nonconventional_T$source_obj <- ifelse(grepl("^DIS_", colnames(nonconventional_T)), "disease", "normal")
nonconventional_T_list <- SplitObject(nonconventional_T, split.by = "source_obj")
names(nonconventional_T_list)

nonconventional_T_list <- lapply(nonconventional_T_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) #<<<
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000) #<<<
  return(x)
})

features <- SelectIntegrationFeatures(object.list = nonconventional_T_list, nfeatures = 3000) #<<<

nonconventional_T_list <- lapply(nonconventional_T_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

anchors <- FindIntegrationAnchors(
  object.list = nonconventional_T_list,
  anchor.features = features,
  reduction = "rpca", #<<<
  dims = 1:30 #<<<
)

nonconventional_T <- IntegrateData(
  anchorset = anchors,
  dims = 1:30 #<<<
)

nonconventional_T <- ScaleData(nonconventional_T, verbose = FALSE)
nonconventional_T <- RunPCA(nonconventional_T, npcs = 50, verbose = TRUE)
ElbowPlot(nonconventional_T) #3x4 #20260324_B_subcluster_PCA_elbowplot
nonconventional_T <- FindNeighbors(nonconventional_T, dims = 1:10) #<<< 
nonconventional_T <- FindClusters(nonconventional_T, resolution = 0.2) #<<<
nonconventional_T <- RunUMAP(nonconventional_T, dims = 1:10) #<<<
DimPlot(nonconventional_T, reduction='umap',label=T)#3x4
DimPlot(nonconventional_T, reduction='umap',split.by='cohort_l2',label=T)
DimPlot(subset(nonconventional_T, subset = manual.cluster == 'naive/TCM CD4+T'), reduction='umap', group.by = 'predicted.celltype.l2',label=T)#3x4
markers <- FindAllMarkers(nonconventional_T, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

nonconventional_T$manual.cluster_l2 <- ifelse(nonconventional_T$seurat_clusters==0,'residual_conventional_naive_T',
                                        ifelse(nonconventional_T$seurat_clusters %in% c(1,2), 'gamma_delta_T',
                                               ifelse(nonconventional_T$seurat_clusters==3,'MAIT',
                                                      ifelse(nonconventional_T$seurat_clusters==4,'immediate early response',
                                                            'contam'))))

nonconventional_T$manual.cluster_l3 <- ifelse(nonconventional_T$seurat_clusters==0,'residual_conventional_naive_T',
                                              ifelse(nonconventional_T$seurat_clusters ==1, 'TRDV1 like gamma_delta_T',
                                                     ifelse(nonconventional_T$seurat_clusters ==2,'TRDV2 TRGV9 gamma_delta_T',
                                                     ifelse(nonconventional_T$seurat_clusters==3,'MAIT',
                                                            ifelse(nonconventional_T$seurat_clusters==4,'immediate early response',
                                                                   'contam')))))

DimPlot(nonconventional_T, group.by='manual.cluster_l3',label=T)#4x6.5 #20260401_nonconventional_T_umap
saveRDS(nonconventional_T, file='data/20260309_pilot/results/version2/nonconventional_T/nonconventional_T_subset.rds')

cell_count_long <- nonconventional_T@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::filter(cohort_l2 == "Lee_p1_base") %>%
  count(sample, manual.cluster_l3, name = "n_cells") %>%
  arrange(sample, manual.cluster_l3)

cell_count_long
nrow(cell_count_long[cell_count_long$manual.cluster_l3=='TRDV1 like gamma_delta_T' ,])
nrow(cell_count_long[cell_count_long$manual.cluster_l3=='TRDV1 like gamma_delta_T' & cell_count_long$n_cells >= 10,]) #78
nrow(cell_count_long[cell_count_long$manual.cluster_l3=='TRDV2 TRGV9 gamma_delta_T' ,])
nrow(cell_count_long[cell_count_long$manual.cluster_l3=='TRDV2 TRGV9 gamma_delta_T' & cell_count_long$n_cells >= 10,]) #51
nrow(cell_count_long[cell_count_long$manual.cluster_l3=='TRDV2 TRGV9 gamma_delta_T' & cell_count_long$n_cells >= 20,])

ggplot(cell_count_long, aes(x=n_cells,fill=manual.cluster_l3))+geom_histogram()+facet_wrap(.~manual.cluster_l3, scales='free')
##
# centroid based continuum - TRDV1, TRDV2/TRGV9
# ilr x
# cellrank - TRDV1, TRDV2/TRGV9
# pseudobulk based -curated/de novo - gdT
##
#####
# =========================================================
# 0. USER CONFIG
# =========================================================
TARGET_ANALYSIS <- "nonconventional_T"   # "Monocyte" or "NK"

# optional: reload clinical metadata from original meta file
ADD_CLINICAL_META <- TRUE
META_PATH <- "data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt"

# common metadata columns in your Seurat object
SAMPLE_COL        <- "sample"
COHORT_COL        <- "cohort_l2"
RESPONSE_COL      <- "binarized_response"
ASSAY_USE         <- "RNA"
TARGET_COHORT     <- "Lee_p1_base"

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
SUBTYPE_COL       <- "manual.cluster_l3"
gdT <- subset(nonconventional_T, subset = manual.cluster_l2 == 'gamma_delta_T')
gdT <- subset(gdT, subset = cohort_l2 == 'Lee_p1_base')

DefaultAssay(gdT) <- 'RNA'
gdT<- JoinLayers(gdT)

gdT <- NormalizeData(gdT, normalization.method = "LogNormalize", scale.factor = 10000) 
gdT <- FindVariableFeatures(gdT, selection.method = "vst", nfeatures = 3000) #<<<

gdT <- ScaleData(gdT, verbose = FALSE)
gdT <- RunPCA(gdT, npcs = 50, verbose = TRUE)

DimPlot(gdT, reduction='pca', group.by='manual.cluster_l3', split.by='')
ElbowPlot(gdT) #3x4 #20260330_CD8T_naive_cytotoxic_PCA_elbow
DimHeatmap(gdT, dims = 1:10, cells = 500, balanced = TRUE) #8x8 #20260401_gdT_PCA_dimheat
# PC: 3,4
ANCHOR_STATE_LOW  <- 'TRDV2 TRGV9 gamma_delta_T' 
ANCHOR_STATE_HIGH <- 'TRDV1 like gamma_delta_T'
CENTROID_PREFIX   <-"gdT_circulating_to_NKlike"

centroid_outdir <- file.path(OUTROOT, "centroid_from_clean_subset")
dir.create(centroid_outdir, recursive = TRUE, showWarnings = FALSE)

continuum_df <- NULL
if (!is.na(ANCHOR_STATE_LOW) && !is.na(ANCHOR_STATE_HIGH)) {
  continuum_res <- calc_two_anchor_score(
    seu = gdT,
    sample_col = "sample",
    celltype_col = SUBTYPE_COL,
    cohort_col = "cohort_l2",
    target_cohort = "Lee_p1_base",
    low_label = ANCHOR_STATE_LOW,
    high_label = ANCHOR_STATE_HIGH,
    assay_use = "RNA",
    rerun_pca = TRUE,
    dims_use = c(3,4),
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
# 3. cellrank
########


if (DO_CELLRANK_EXPORT) {
  cellrank_outdir <- file.path(OUTROOT, "CellRank_input_from_clean_subset")
  export_cellrank_input(
    seu = gdT,
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

trdv2_trgv9_marker <- c('TRDV2','TRGV9','TRDC')#2
trdv1_marker <- c('TRDV1','KLRC2', 'KLRF1','CD8B','TIGIT','CD8A','IKZF2','FCRL6')#1
trdv1.markers <- FindMarkers(gdT, ident.1 = 1, ident.2 = 2)
trdv1.markers <- trdv1.markers[trdv1.markers$avg_log2FC >0,]
head(trdv1.markers[c(order(trdv1.markers$p_val_adj)),],20)

VlnPlot(gdT, features=trdv1_marker,pt.size=0)
VlnPlot(gdT, features=trdv2_trgv9_marker,pt.size=0) # cluster 1= resting

#########
# 4. curated gene set
########
SUBTYPE_COL       <- "manual.cluster_l2"
get_analysis_defaults <- function(target_analysis) {
  norm_symbol <- function(x) {
    unique(toupper(trimws(as.character(x))))
  }
  
  if (target_analysis == "nonconventional_T") {
    
    list(
      gene_blocklist_patterns = c(
        "^MT-", "^RPS", "^RPL", "^HSP",
        "^IG[HKL]",
        "^TRA", "^TRB", "^TRD", "^TRG"
      ),
      custom_gene_sets = list(
        GDT_MEMORY_LIKE = c(
          "IL7R", "CCR7", "TCF7", "LEF1", "SELL",
          "LTB", "MAL", "RCAN3", "AQP3", "GPR183",
          "INPP4B", "CD27", "TMEM123", "PRKCQ-AS1"
        ),
        
        GDT_CYTOTOXIC_EFFECTOR = c(
          "NKG7", "PRF1", "GZMB", "GZMH", "GNLY",
          "CCL4", "CCL5", "CTSW", "FGFBP2", "KLRD1",
          "KLRK1", "PLEK", "ADGRG1"
        ),
        
        GDT_TERMINAL_NKLIKE = c(
          "FCGR3A", "ADGRG1", "KLRC2", "KLRC3", "KLRF1",
          "KIR2DL3", "FCRL6", "ZEB2", "S1PR5", "AOAH",
          "CMC1", "SPON2", "FGFBP2"
        ),
        
        GDT_ADAPTIVE_VD1_LIKE = c(
          "KLRC2", "KLRC3", "TIGIT", "FCRL6", "CD8A",
          "CD8B", "ZEB2", "AOAH", "S1PR5", "ADGRG1",
          "LGALS1", "IFITM2", "HLA-DRB1"
        ),
        
        GDT_BLOOD_VG9VD2_LIKE = c(
          "KLRG1", "FGFBP2", "GNLY", "CCL4", "GZMB",
          "PRF1", "SPON2", "PTGDS", "KLRD1", "KLRK1",
          "KLRC1", "PLEK"
        ),
        
        GDT_TCR_ACTIVATION_IMMEDIATE = c(
          "CD69", "NR4A1", "NR4A2", "NR4A3",
          "FOS", "FOSB", "JUN", "JUNB",
          "EGR1", "EGR2", "TNFAIP3", "NFKBIA",
          "DUSP1", "ZFP36", "PPP1R15A", "CSRNP1"
        ),
        
        GDT_IFN_RESPONSE = c(
          "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6",
          "IFI44", "IFI44L", "MX1", "MX2", "OAS1",
          "OASL", "STAT1", "IRF7", "GBP1", "IFITM1", "IFITM2"
        ),
        
        GDT_CHECKPOINT_INHIBITORY = c(
          "TIGIT", "PDCD1", "LAG3", "HAVCR2", "CTLA4",
          "IKZF2", "LGALS1", "CD160", "TOX", "TOX2",
          "ENTPD1", "BATF", "PRDM1"
        ),
        
        GDT_HLAII_ACTIVATED = c(
          "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1",
          "HLA-DQB1", "CD74", "CIITA", "GZMK", "TIGIT"
        ),
        
        GDT_PROLIFERATION_CYCLING = c(
          "MKI67", "TYMS", "PCLAF", "RRM2", "TK1",
          "STMN1", "HMGB2", "BIRC5", "MCM5", "MCM6",
          "PTTG1", "TOP2A"
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
        "GOBP_REGULATION_OF_T_CELL_ACTIVATION",
        "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
        "GOBP_LYMPHOCYTE_MIGRATION",
        "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
        "GOBP_POSITIVE_REGULATION_OF_CELL_KILLING",
        "GOBP_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY",
        "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
        "GOBP_CYTOKINE_PRODUCTION",
        "REACTOME_TCR_SIGNALING",
        "REACTOME_DOWNSTREAM_TCR_SIGNALING",
        "REACTOME_SIGNALING_BY_INTERLEUKINS",
        "REACTOME_PD_1_SIGNALING",
        "REACTOME_INTERFERON_SIGNALING"
      )
    )
    
  } else {
    stop("TARGET_ANALYSIS must be 'CD8T' or 'gdT'.")
  }
}
cfg <- get_analysis_defaults(TARGET_ANALYSIS)
GENE_BLOCKLIST_PATTERNS <- cfg$gene_blocklist_patterns
CUSTOM_GENE_SETS <- cfg$custom_gene_sets
MSIGDB_TARGETS <- cfg$msigdb_targets

curated_outdir <- file.path(OUTROOT, "curated_program_scores_from_clean_subset")
curated_res <- run_curated_pseudobulk_scores(
  seu_use = gdT,
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
k_vec <- 3:5
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
meta1 <- gdT@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

keep_cells_non_na <- meta1$cell[!is.na(meta1[[sample_col]]) & meta1[[sample_col]] != ""]
seu_non_na <- subset(gdT, cells = keep_cells_non_na)

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
ph #8x9 #

nmf_outdir <- file.path(OUTROOT, "geneNMF_from_clean_subset")
nmf_res <- run_gene_nmf_features(
  seu_use = gdT,
  outdir = nmf_outdir,
  sample_col = SAMPLE_COL,
  cohort_col = COHORT_COL,
  target_cohort = TARGET_COHORT,
  assay_use = ASSAY_USE,
  min_cells_per_sample = MIN_CELLS_PB,
  k_vec = 3:5,
  nfeatures_hvg = 2000,
  nmf_seed = 123,
  gene_blocklist_patterns = GENE_BLOCKLIST_PATTERNS
)
geneNMF_metaprograms <- readRDS('data/20260309_pilot/results/version2/nonconventional_T/geneNMF_from_clean_subset/05_geneNMF_metaprograms.rds')
mp_genes <- geneNMF_metaprograms$metaprograms.genes
mp_metrics <- geneNMF_metaprograms$metaprograms.metrics

anno_colors <- brewer.pal(n=5, name="Paired")
names(anno_colors) <- names(geneNMF_metaprograms$metaprograms.genes)
ph <- plotMetaPrograms(geneNMF_metaprograms, annotation_colors = anno_colors)
ph #20260401_gdT_NMF_metaprogam_heatmap

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

CELLRANK_PATIENT_CSV = 'data/20260309_pilot/results/version2/nonconventional_T/CellRank_output_from_clean_subset/cellrank_b3_patient_features.csv'
if (nzchar(CELLRANK_PATIENT_CSV) && file.exists(CELLRANK_PATIENT_CSV)) {
  final_feature_df_cellrank <- merge_cellrank_features(final_feature_df, CELLRANK_PATIENT_CSV)
  write.csv(final_feature_df_cellrank, file.path(OUTROOT, paste0("final_feature_df_with_cellrank_", sanitize_name(TARGET_ANALYSIS), ".csv")), row.names = FALSE)
  saveRDS(final_feature_df_cellrank, file.path(OUTROOT, paste0("final_feature_df_with_cellrank_", sanitize_name(TARGET_ANALYSIS), ".rds")))
}

final_feature_df_cellrank_filt <- final_feature_df_cellrank[,!(colnames(final_feature_df_cellrank) %in% c('denovo_gene_mp1__singscore','denovo_gene_mp3__singscore','denovo_gene_mp4__singscore','denovo_gene_mp1__n_overlap','denovo_gene_mp2__n_overlap', 'denovo_gene_mp3__n_overlap', 'denovo_gene_mp4__n_overlap', 'Centroid_n_score_cells',  'Centroid_keep_for_analysis', 'cellrank_cr_b3_priming_mean', 'cellrank_cr_b3_priming_q75','cellrank_qc_b3_subset_frac_TRDV1 like gamma_delta_T','cellrank_qc_b3_subset_frac_TRDV2 TRGV9 gamma_delta_T'))]
OUT_CSV <- "data/20260309_pilot/results/version2/nonconventional_T/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv"
fwrite(as.data.table(final_feature_df_cellrank_filt), OUT_CSV)

