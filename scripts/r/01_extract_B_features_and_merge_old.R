# 2026.03.24
# make feature library in B

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

## cell type proportion
# n_cells per patients celltype >= 10

## pesudobulk level
# n_cells per patients celltype >= 20

## 
cell_df <- b_sbust@meta.data[b_sbust@meta.data$time =='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  dplyr::select(sample, manual.cluster_l2)
count_df <- b_sbust@meta.data[b_sbust@meta.data$time =='base' & b_sbust@meta.data$cohort =='Lee_p1',] %>%
  count(sample, manual.cluster_l2, name = "n_cells")
total_df <- b_sbust@meta.data[b_sbust@meta.data$time =='base' & b_sbust@meta.data$cohort =='Lee_p1',] %>%
  count(sample, name = "total_cells")
prop_df <- count_df %>%
  left_join(total_df, by = "sample") %>%
  mutate(prop = n_cells / total_cells)
ggplot(prop_df[prop_df$manual.cluster_l2=='intermediate_resting',], aes(x=1, y=n_cells))+geom_boxplot()+geom_beeswarm()+ylim(c(0,25))
ggplot(prop_df, aes(x='celltype', y=n_cells))+geom_boxplot(outlier.size=0)+geom_beeswarm(size =1, alpha=0.6,cex=4)+
  facet_wrap(.~manual.cluster_l2,scales='free_y') #20260324_b_subset_patient_n_cells #5x6
ggplot(prop_df[prop_df$sample %in% selected_patients$V1,], aes(x='celltype', y=n_cells))+geom_boxplot(outlier.size=0)+geom_beeswarm(size =1, alpha=0.6,cex=4)+
  facet_wrap(.~manual.cluster_l2,scales='free_y')
for (i in unique(prop_df$manual.cluster)){
  print(paste0(i,' ',min(prop_df[prop_df$sample %in% selected_patients$V1 & prop_df$manual.cluster==i,]$n_cells)))
}
# [1] "naive/TCM CD4+T 135"
# [1] "CD14 Mono 100"
# [1] "effector CD8+T 60"
# [1] "NK 31"
# [1] "CD16 Mono 11"
# [1] "T/NK 12"
# [1] "B 11"
# [1] "ambiguous 1"
# [1] "cDC2 1"
# [1] "navie CD8+T 1"
# [1] "platelet 2"

## cell type annotation
DimPlot(seu_integrated, reduction = "umap", group.by='manual.cluster',split.by='manual.cluster')+facet_wrap(~manual.cluster)
DimPlot(seu_integrated, reduction = "umap", group.by='integrated_snn_res.0.3')
DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l2',label=T)

seu_integrated$plot_group <- ifelse(seu_integrated$manual.cluster == 'NK','target','other')
seu_integrated$plot_group <- factor(seu_integrated$plot_group, levels=c('other','target'))
DimPlot(
  seu_integrated, reduction='umap',
  group.by = 'plot_group',cols=c('other'='lightgrey','target'='navy'),
  order='target'
)#3x4 #20260324_cite_gpT_dimplot #20260324_cite_Treg_dimplot #20260324_cite_MAIT_dimplot #20260324_cite_NK_CD56bright_dimplot #20260324_cite_NK_CD8TEM_dimplot #20260324_cite_CD8TCM_dimplot #20260324_cite_CD8Naive_dimplot #20260324_cite_CD4Naive_dimplot #20260324_cite_CD4TCM_dimplot #20260324_cite_CD4TEM_dimplot #20260324_cite_CD4CTL_dimplot #20260324_cite_NK_dimplot #20260324_cite_B_naive_dimplot #20260324_cite_B_intermediate_dimplot #20260324_cite_B_memory_dimplot #20260324_manual_Naive_TCM_CD4T_dimplot #20260324_manual_Naive_CD8T_dimplot #20260324_manual_effector_CD8T_dimplot #20260324_manual_T_NK_dimplot #20260324_manual_NK_dimplot

##############
## B cell subclustering

B_subset<- subset(seu_integrated, idents = "B")
DefaultAssay(B_subset) <- 'RNA'
B_subset$source_obj <- ifelse(grepl("^DIS_", colnames(B_subset)), "disease", "normal")
b_list <- SplitObject(B_subset, split.by = "source_obj")
names(b_list)

b_list <- lapply(b_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) #<<<
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000) #<<<
  return(x)
})

features <- SelectIntegrationFeatures(object.list = b_list, nfeatures = 3000) #<<<

b_list <- lapply(b_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

anchors <- FindIntegrationAnchors(
  object.list = b_list,
  anchor.features = features,
  reduction = "rpca", #<<<
  dims = 1:30 #<<<
)

b_sbust <- IntegrateData(
  anchorset = anchors,
  dims = 1:30 #<<<
)

b_sbust <- ScaleData(b_sbust, verbose = FALSE)
b_sbust <- RunPCA(b_sbust, npcs = 50, verbose = TRUE)
DimPlot(b_sbust, reduction='pca', group.by='predicted.celltype.l2')
DimPlot(b_sbust, reduction='pca', group.by='manual.cluster_l2')
ElbowPlot(b_sbust) #3x4 #20260324_B_subcluster_PCA_elbowplot
DimHeatmap(b_sbust, dims = 1:12, cells = 500, balanced = TRUE) #8x8 #20260324_B_subcluster_PCA_pcaheatmap
b_sbust <- FindNeighbors(b_sbust, dims = 1:15) #<<< 
b_sbust <- FindClusters(b_sbust, resolution = 0.1) #<<<
b_sbust <- RunUMAP(b_sbust, dims = 1:15) #<<<
DimPlot(b_sbust, reduction='umap')#3x4
DimPlot(subset(b_sbust, predicted.celltype.l2 %in% c('B intermediate','B memory','B naive')), reduction='umap',group.by='predicted.celltype.l2', split.by = 'predicted.celltype.l2')#4x8 #20260324_b_subset_cite_l2
markers <- FindAllMarkers(b_sbust, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(b_sbust, features = top10$gene) + NoLegend()
# navie: TCL1A, IGHD, FCER2, IL4R
# memory (activated/atypical/inflamed) - GPR183, CLECL1, HMOX1, SOX5
# plasma: XBP1, FKBP11, TXNDC5, PDIA4, LMAN1, MANF
b_sbust$manual.cluster_l2 <- ifelse(b_sbust$seurat_clusters==0,'naive',
                                    ifelse(b_sbust$seurat_clusters==4, 'intermediate_resting',
                                           ifelse(b_sbust$seurat_clusters==1,'memory',
                                                  ifelse(b_sbust$seurat_clusters==5,'plasma','contam'))))
FeaturePlot(b_sbust, features = c('CD27','TNFRSF13B','CD80','CD86','BANK1'), order =T)
ggplot(b_sbust@meta.data, aes(x=manual.cluster_l2, fill=cohort_l2))+geom_bar(position='dodge')
DimPlot(b_sbust, reduction='umap', group.by='manual.cluster_l2',label=T)#6x4 #20260324_B_subcluster
saveRDS(b_sbust, file='data/20260309_pilot/results/version2/B/b_sbust.rds')

### remove contamination, plasma, intermediate_resting
bad_clusters <- c("plasma", "contam", "intermediate_resting")
B_clean <- subset(b_sbust, subset = !(manual.cluster_l2 %in% bad_clusters))
DefaultAssay(B_clean) <- 'RNA'
b_list <- SplitObject(B_clean, split.by = "source_obj")
b_list <- lapply(b_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  x
})
features <- SelectIntegrationFeatures(object.list = b_list, nfeatures = 3000)
b_list <- lapply(b_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  x
})
anchors <- FindIntegrationAnchors(
  object.list = b_list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30
)
b_subust_clean <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)
b_subust_clean <- ScaleData(b_subust_clean, verbose = FALSE)
b_subust_clean <- RunPCA(b_subust_clean, npcs = 50, verbose = TRUE)
ElbowPlot(b_subust_clean) #3x4 #20260324_B_subcluster_PCA_elbowplot
DimHeatmap(b_subust_clean, dims = 1:12, cells = 500, balanced = TRUE) #8x8 #20260324_B_subcluster_PCA_pcaheatmap

################################
####### patient level propotion
celltypes_use <- c("intermediate_resting", "memory", "naive")
cohorts_use   <- c("Lee_p1_base_DCB","Lee_p1_base_NCB", "Lee_p1_1st_DCB","Lee_p1_1st_NCB", "AIDA__NA")
meta_df <- b_sbust@meta.data %>%
  as.data.frame()
meta_df$sample <- ifelse(is.na(meta_df$sample), meta_df$donor_id,meta_df$sample)
meta_df$cohort_l2_resp <- paste0(meta_df$cohort_l2,'_',meta_df$binarized_response)

plot_df <- meta_df %>%
  filter(
    cohort_l2_resp %in% cohorts_use,
    manual.cluster_l2 %in% celltypes_use
  ) %>%
  count(cohort_l2_resp, sample, manual.cluster_l2, name = "n") %>%
  group_by(cohort_l2_resp, sample) %>%
  complete(manual.cluster_l2 = celltypes_use, fill = list(n = 0)) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

naive_df <- plot_df %>%
  filter(manual.cluster_l2 == "naive") %>%
  select(cohort_l2_resp, sample, naive_prop = prop)

base_order_df <- naive_df %>%
  filter(cohort_l2_resp %in% c("Lee_p1_base_DCB","Lee_p1_base_NCB")) %>%
  arrange(naive_prop) %>%
  mutate(base_rank = row_number()) %>%
  select(sample, base_rank)

# - base: 자기 순서
# - 1st : base 순서 따라감
# - AIDA: AIDA 자체 순서
plot_order_df <- naive_df %>%
  left_join(base_order_df, by = "sample") %>%
  group_by(cohort_l2_resp) %>%
  arrange(naive_prop, .by_group = TRUE) %>%
  mutate(self_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    plot_rank = case_when(
      cohort_l2_resp %in% c("Lee_p1_base_DCB","Lee_p1_base_NCB") ~ base_rank,
      cohort_l2_resp %in% c("Lee_p1_1st_DCB","Lee_p1_1st_NCB")  ~ base_rank,
      cohort_l2_resp == "AIDA__NA"        ~ self_rank
    )
  )

# 5) x축용 cohort별 unique label 생성
plot_df2 <- plot_df %>%
  left_join(plot_order_df %>% select(cohort_l2_resp, sample, plot_rank),
            by = c("cohort_l2_resp", "sample")) %>%
  mutate(sample_plot = paste(cohort_l2_resp, sample, sep = "___"))

# facet 순서까지 고정하고 싶으면
cohort_levels <- c("Lee_p1_base_DCB", "Lee_p1_1st_DCB","Lee_p1_base_NCB","Lee_p1_1st_NCB", "AIDA__NA")
plot_df2$cohort_l2_resp <- factor(plot_df2$cohort_l2_resp, levels = cohort_levels)

sample_levels <- plot_df2 %>%
  distinct(cohort_l2_resp, sample, sample_plot, plot_rank) %>%
  arrange(cohort_l2_resp, plot_rank, sample) %>%
  pull(sample_plot)

plot_df2$sample_plot <- factor(plot_df2$sample_plot, levels = sample_levels)
ggplot(plot_df2, aes(x = sample_plot, y = prop, fill = factor(manual.cluster_l2,levels=c('memory','intermediate_resting','naive')))) +
  geom_col(width = 0.9) +
  facet_grid(. ~ factor(cohort_l2_resp,levels=c('AIDA__NA',"Lee_p1_base_DCB", "Lee_p1_1st_DCB","Lee_p1_base_NCB","Lee_p1_1st_NCB")), scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c(
    "intermediate_resting" = "#E69F00",
    "memory" = "#56B4E9",
    "naive" = "navy"
  )) +
  labs(
    x = "Sample",
    y = "Proportion within B cells",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=2),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )#4x10 #20260324_B_subcluster_change


################################


############################
## A. composition/balance level
############################
# 1. centroid (anchor) based scoring
# #
library(Seurat)
library(dplyr)
library(tibble)

calc_b_naive_memory_score <- function(
    b_subust_clean,
    sample_col = "sample",
    celltype_col = "cluster_l2",
    naive_label = "naive",
    memory_label = "memory",
    assay_use = "RNA",
    rerun_pca = TRUE,
    reduction_use = "pca",
    dims_use = 1:4,
    npcs = 20,
    min_cells_per_patient = 10,
    score_cut = 0.5,
    eps = 1e-8
) {
  stopifnot(inherits(b_subust_clean, "Seurat"))
  stopifnot(sample_col %in% colnames(b_subust_clean@meta.data))
  stopifnot(celltype_col %in% colnames(b_subust_clean@meta.data))
  
  meta_clean <- b_subust_clean@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    mutate(
      .sample = .data[[sample_col]],
      .celltype = .data[[celltype_col]]
    )
  
  keep_ok <- unique(meta_clean$.celltype)
  if (!all(keep_ok %in% c(naive_label, memory_label))) {
    stop("b_subust_clean should contain only naive_label and memory_label cells.")
  }
  
  if (sum(meta_clean$.celltype == naive_label) < 2) {
    stop("Too few naive anchor cells in b_subust_clean.")
  }
  if (sum(meta_clean$.celltype == memory_label) < 2) {
    stop("Too few memory anchor cells in b_subust_clean.")
  }
  
  seu_work <- b_subust_clean
  DefaultAssay(seu_work) <- assay_use
  
  if (rerun_pca) {
    seu_work <- NormalizeData(seu_work, verbose = FALSE)
    seu_work <- FindVariableFeatures(seu_work, verbose = FALSE)
    seu_work <- ScaleData(seu_work, verbose = FALSE)
    seu_work <- RunPCA(seu_work, npcs = npcs, verbose = FALSE)
    reduction_use <- "pca"
  } else {
    if (!reduction_use %in% Reductions(seu_work)) {
      stop(paste0("Reduction '", reduction_use, "' not found in b_subust_clean."))
    }
  }
  
  emb <- Embeddings(seu_work, reduction = reduction_use)
  if (max(dims_use) > ncol(emb)) {
    stop("dims_use exceeds available dimensions.")
  }
  emb <- emb[, dims_use, drop = FALSE]
  
  df_clean <- emb %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    inner_join(meta_clean, by = "cell")
  
  dim_cols <- colnames(emb)
  
  naive_centroid <- colMeans(df_clean[df_clean$.celltype == naive_label, dim_cols, drop = FALSE])
  memory_centroid <- colMeans(df_clean[df_clean$.celltype == memory_label, dim_cols, drop = FALSE])
  
  v <- memory_centroid - naive_centroid
  vv <- sum(v^2)
  if (vv == 0) {
    stop("Naive and memory centroids are identical in the selected latent space.")
  }
  
  X <- as.matrix(df_clean[, dim_cols, drop = FALSE])
  
  raw_score <- as.numeric(
    (X - matrix(naive_centroid, nrow = nrow(X), ncol = length(naive_centroid), byrow = TRUE)) %*% v / vv
  )
  
  clipped_score <- pmin(pmax(raw_score, 0), 1)
  
  df_clean$naive_memory_score_raw <- raw_score
  df_clean$naive_memory_score <- clipped_score
  df_clean$naive_memory_bin <- ifelse(df_clean$naive_memory_score < score_cut, "naive_like", "memory_like")
  
  nm_summary <- df_clean %>%
    group_by(.sample) %>%
    summarise(
      n_score_cells = n(),
      mean_score = mean(naive_memory_score),
      median_score = median(naive_memory_score),
      sd_score = sd(naive_memory_score),
      iqr_score = IQR(naive_memory_score),
      
      frac_naive_like = mean(naive_memory_bin == "naive_like"),
      frac_memory_like = mean(naive_memory_bin == "memory_like"),
      
      observed_frac_naive = mean(.celltype == naive_label),
      observed_frac_memory = mean(.celltype == memory_label),
      
      .groups = "drop"
    ) %>%
    mutate(
      log2ratio_memory_naive = log2((frac_memory_like + eps) / (frac_naive_like + eps)),
      skew_memory_minus_naive = frac_memory_like - frac_naive_like,
      keep_for_analysis = n_score_cells >= min_cells_per_patient
    ) %>%
    arrange(median_score)
  
  clean_meta_add <- df_clean %>%
    select(
      cell,
      B_naive_memory_score = naive_memory_score,
      B_naive_memory_score_raw = naive_memory_score_raw,
      B_naive_memory_bin = naive_memory_bin
    )
  
  b_subset_clean_scored <- b_subust_clean
  clean_meta_add2 <- clean_meta_add[match(Cells(b_subset_clean_scored), clean_meta_add$cell), -1, drop = FALSE]
  rownames(clean_meta_add2) <- Cells(b_subset_clean_scored)
  b_subset_clean_scored <- AddMetaData(b_subset_clean_scored, metadata = clean_meta_add2)
  
  return(list(
    b_subset_clean_scored = b_subset_clean_scored,
    cell_scores = df_clean,
    patient_summary = nm_summary,
    naive_centroid = naive_centroid,
    memory_centroid = memory_centroid,
    score_cut = score_cut,
    reduction_used = reduction_use,
    dims_used = dims_use
  ))
}
b_subust_clean_poi <- subset(b_subust_clean, cohort_l2 == 'Lee_p1_base')
res <- calc_b_continuum_score(
  seu = b_subust_clean_poi,
  sample_col = "sample",
  celltype_col = "manual.cluster_l2",
  naive_label = "naive",
  intermediate_label = "intermediate_resting",
  memory_label = "memory",
  assay_use = "RNA",      # 필요시 "SCT"로 변경
  rerun_pca = TRUE,
  dims_use = 1:4,
  npcs = 20,
  min_cells_per_patient = 10
)

head(res$patient_summary)
patient_summary_melt <- reshape2::melt(res$patient_summary, id.vars=c('.sample','n_b_cells'),measure.vars=c('mean_score','median_score','sd_score','iqr_score','frac_naive_like','frac_memory_like','observed_frac_naive','observed_frac_memory','skew_memory_minus_naive','log2ratio_memory_naive','entropy_norm'))
ggplot(patient_summary_melt, aes(x=value))+geom_histogram(bins=30)+facet_wrap(.~variable, scales='free')+theme_bw()#20260324_B_naive_memory_continuum_score
table(res$seu_scored$B_continuum_bin, useNA = "ifany")
res$thresholds

saveRDS(res, file='data/20260309_pilot/results/version2/B/naive_memory_continuum_data.rds')
############################
# 2. pseudotime like lineage score


############################
## B. pseudobulk/module level
############################
# 1. curated gene set
# 1.1 inflammation 0 MSigDB Hallmark
# HALLMARK_INTERFERON_ALPHA_RESPONSE,HALLMARK_INTERFERON_GAMMA_RESPONSE, HALLMARK_TNFA_SIGNALING_VIA_NFKB, HALLMARK_IL6_JAK_STAT3_SIGNALING, HALLMARK_INFLAMMATORY_RESPONSE
#CUSTOM_B_AP_MHCII_CORE = c("CD74","HLA-DMA", "HLA-DMB","HLA-DOA", "HLA-DOB","HLA-DPA1", "HLA-DPB1","HLA-DQA1", "HLA-DQA2","HLA-DQB1", "HLA-DQB2","HLA-DRA","HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "CIITA","IFI30", "LGMN","CTSB", "CTSD", "CTSL", "CTSS")
#CUSTOM_B_AP_IFN_INDUCED = c("TAP1", "TAP2", "B2M", "PSMB8", "PSMB9", "HLA-A", "HLA-B", "HLA-C")
#CUSTOM_B_AP_COSTIM = c("CD80", "CD86", "CD40", "ICAM1")

# sample filter: Lee_p1_base, min_cells_per_Sample
# gene filter: min_samples_passing_detection >=2, non-zero cell>=10, pseudobulk logCPM varaince >0
## =========================================================
## 0. packages
## =========================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(edgeR)
library(msigdbr)
library(singscore)
library(ssGSEA2)

## =========================================================
## 1. input 설정
## =========================================================
#seu_use <- b_subust_clean
sample_col <- "sample"
cohort_col <- "cohort_l2"
target_cohort <- "Lee_p1_base"
assay_use <- "RNA"
species_use <- "Homo sapiens"

outdir <- "data/20260309_pilot/results/version2/B/lineage_program_scores_stepwise"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_sample <- 20
min_geneset_overlap <- 3

## =========================================================
## 2. metadata / cohort 상태 확인
## =========================================================
meta0 <- b_subust_clean@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

cat("sample_col 존재 여부:\n")
print(sample_col %in% colnames(meta0))

cat("\ncohort_col 존재 여부:\n")
print(cohort_col %in% colnames(meta0))

cat("\ncohort 분포:\n")
print(table(meta0[[cohort_col]], useNA = "ifany"))

write.csv(
  data.frame(cohort = meta0[[cohort_col]]) %>% count(cohort, name = "n_cells"),
  file.path(outdir, "00_cohort_distribution_raw.csv"),
  row.names = FALSE
)

## =========================================================
## 3. target cohort만 남기기
## =========================================================
keep_cells_cohort <- meta0$cell[!is.na(meta0[[cohort_col]]) & meta0[[cohort_col]] == target_cohort]
seu_cohort <- subset(b_subust_clean, cells = keep_cells_cohort)

cat("\nTarget cohort =", target_cohort, " 에 해당하는 cell 수:\n")
print(ncol(seu_cohort))

meta1 <- seu_cohort@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

cat("\nTarget cohort subset에서 sample 분포 (NA 포함):\n")
print(table(meta1[[sample_col]], useNA = "ifany"))

write.csv(
  data.frame(sample = meta1[[sample_col]]) %>% count(sample, name = "n_cells"),
  file.path(outdir, "01_sample_distribution_in_target_cohort_raw.csv"),
  row.names = FALSE
)

## =========================================================
## 4. sample NA / 빈값 제거
## =========================================================
cat("\nNA sample 수:\n")
print(sum(is.na(meta1[[sample_col]])))

cat("\n빈 문자열 sample 수:\n")
print(sum(meta1[[sample_col]] == "", na.rm = TRUE))

keep_cells_non_na <- meta1$cell[!is.na(meta1[[sample_col]]) & meta1[[sample_col]] != ""]
seu_non_na <- subset(seu_cohort, cells = keep_cells_non_na)

cat("\nTarget cohort + non-NA sample cell 수:\n")
print(ncol(seu_non_na))

meta2 <- seu_non_na@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

sample_meta2 <- meta2 %>%
  transmute(sample = as.character(.data[[sample_col]])) %>%
  count(sample, name = "n_cells") %>%
  arrange(desc(n_cells))

cat("\nTarget cohort 안에서 non-NA sample 분포:\n")
print(sample_meta2)

write.csv(sample_meta2, file.path(outdir, "02_sample_distribution_non_na.csv"), row.names = FALSE)

## =========================================================
## 5. min_cells_per_sample 적용
## =========================================================
keep_samples <- sample_meta2 %>%
  filter(n_cells >= min_cells_per_sample) %>%
  pull(sample)

cat("\nmin_cells_per_sample =", min_cells_per_sample, "조건 통과 sample 수:\n")
print(length(keep_samples)) #69 patients

cat("\n통과 sample 목록:\n")
print(keep_samples)

write.csv(
  data.frame(sample = keep_samples),
  file.path(outdir, "03_keep_samples.csv"),
  row.names = FALSE
)

if (length(keep_samples) < 2) {
  stop("2개 미만 sample만 남았습니다. cohort/sample/min_cells_per_sample를 먼저 점검하세요.")
}

keep_cells <- meta2$cell[meta2[[sample_col]] %in% keep_samples]
seu_pb <- subset(seu_non_na, cells = keep_cells)

cat("\nPseudobulk용 object cell 수:\n")
print(ncol(seu_pb))

## =========================================================
## 5.5 gene filtering (cell-level detection + pseudobulk variance)
## =========================================================
min_nonzero_cells_per_sample <- 10
min_samples_passing_detection <- 2   # 더 보수적으로 하려면 2로 변경

DefaultAssay(seu_pb) <- assay_use

## counts matrix 가져오기
counts_mat <- GetAssayData(seu_pb, assay = assay_use, layer = "counts")

cat("\ncounts matrix dim (genes x cells):\n")
print(dim(counts_mat))

## cell -> sample 매핑
cell_sample_map <- seu_pb@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(.data[[sample_col]])
  )

## counts matrix column 순서와 맞추기
cell_sample_vec <- cell_sample_map$sample[match(colnames(counts_mat), cell_sample_map$cell)]

if (any(is.na(cell_sample_vec))) {
  stop("Some cells in counts_mat could not be mapped to sample labels.")
}

## ---------------------------------------------------------
## 5.5.1 gene별 sample별 non-zero cell 수 계산
## ---------------------------------------------------------
sample_levels <- unique(cell_sample_vec)

nonzero_count_list <- lapply(sample_levels, function(s) {
  idx <- which(cell_sample_vec == s)
  if (length(idx) == 1) {
    nz <- as.numeric(counts_mat[, idx] > 0)
  } else {
    nz <- Matrix::rowSums(counts_mat[, idx, drop = FALSE] > 0)
  }
  data.frame(
    gene = rownames(counts_mat),
    sample = s,
    nonzero_cells = nz,
    stringsAsFactors = FALSE
  )
})

gene_detection_long <- bind_rows(nonzero_count_list)

## sample별 detection matrix (genes x samples)로 변환
gene_detection_wide <- gene_detection_long %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = nonzero_cells
  )

gene_detection_mat <- gene_detection_wide %>%
  column_to_rownames("gene") %>%
  as.matrix()

cat("\ngene detection matrix dim (genes x samples):\n")
print(dim(gene_detection_mat))

## 각 gene이 몇 개 sample에서 min_nonzero_cells_per_sample 이상 검출되는지
n_samples_detected <- rowSums(gene_detection_mat >= min_nonzero_cells_per_sample)

gene_filter_detection_df <- data.frame(
  gene = rownames(gene_detection_mat),
  n_samples_detected = n_samples_detected,
  pass_detection = n_samples_detected >= min_samples_passing_detection,
  stringsAsFactors = FALSE
)

cat("\nDetection filter 통과 gene 수:\n")
print(sum(gene_filter_detection_df$pass_detection)) #10840

write.csv(
  gene_filter_detection_df,
  file.path(outdir, "03a_gene_filter_detection.csv"),
  row.names = FALSE
)

## detection 기준 통과 gene
genes_pass_detection <- gene_filter_detection_df$gene[gene_filter_detection_df$pass_detection]

## ---------------------------------------------------------
## 5.5.2 detection 기준 통과 gene만 사용해서 pseudobulk 생성
## ---------------------------------------------------------
seu_pb$.pb_sample <- as.character(seu_pb@meta.data[[sample_col]])
pb_list <- AggregateExpression(
  object = seu_pb,
  assays = assay_use,
  group.by = ".pb_sample",
  return.seurat = FALSE,
  verbose = FALSE
)

pb_counts_all <- as.matrix(pb_list[[assay_use]])
pb_counts_all <- pb_counts_all[, keep_samples[keep_samples %in% colnames(pb_counts_all)], drop = FALSE]

rownames(pb_counts_all) <- toupper(rownames(pb_counts_all))
pb_counts_all <- rowsum(pb_counts_all, group = rownames(pb_counts_all), reorder = FALSE)
pb_counts_all <- as.matrix(pb_counts_all)

## detection gene도 uppercase 맞춤
genes_pass_detection <- toupper(genes_pass_detection)
genes_pass_detection <- intersect(genes_pass_detection, rownames(pb_counts_all))

pb_counts_det <- pb_counts_all[genes_pass_detection, , drop = FALSE]

cat("\nDetection-filtered pseudobulk dim:\n")
print(dim(pb_counts_det))

## ---------------------------------------------------------
## 5.5.3 pseudobulk variance filter
## ---------------------------------------------------------
pb_logcpm_det <- edgeR::cpm(pb_counts_det, log = TRUE, prior.count = 1)

gene_var <- apply(pb_logcpm_det, 1, var, na.rm = TRUE)

gene_filter_var_df <- data.frame(
  gene = rownames(pb_logcpm_det),
  variance_logcpm = gene_var,
  pass_variance = gene_var > 0,
  stringsAsFactors = FALSE
)

cat("\nVariance > 0 filter 통과 gene 수:\n")
print(sum(gene_filter_var_df$pass_variance))

write.csv(
  gene_filter_var_df,
  file.path(outdir, "03b_gene_filter_variance.csv"),
  row.names = FALSE
)

genes_pass_variance <- gene_filter_var_df$gene[gene_filter_var_df$pass_variance]

## ---------------------------------------------------------
## 5.5.4 최종 gene set
## ---------------------------------------------------------
genes_keep_final <- intersect(genes_pass_detection, genes_pass_variance)

cat("\n최종 filtering 통과 gene 수:\n")
print(length(genes_keep_final))

write.csv(
  data.frame(gene = genes_keep_final),
  file.path(outdir, "03c_gene_filter_final_keep_genes.csv"),
  row.names = FALSE
)

## 최종 pseudobulk counts / logCPM
pb_counts2 <- pb_counts_all[genes_keep_final, , drop = FALSE]
pb_logcpm <- edgeR::cpm(pb_counts2, log = TRUE, prior.count = 1)

cat("\nFinal filtered pseudobulk counts dim:\n")
print(dim(pb_counts2))

cat("\nFinal filtered pseudobulk logCPM dim:\n")
print(dim(pb_logcpm))

saveRDS(pb_counts2, file.path(outdir, "04_pseudobulk_counts_filtered.rds"))
saveRDS(pb_logcpm, file.path(outdir, "05_pseudobulk_logcpm_filtered.rds"))

write.csv(
  pb_counts2[1:min(20, nrow(pb_counts2)), , drop = FALSE],
  file.path(outdir, "04_pseudobulk_counts_filtered_head20.csv")
)

write.csv(
  pb_logcpm[1:min(20, nrow(pb_logcpm)), , drop = FALSE],
  file.path(outdir, "05_pseudobulk_logcpm_filtered_head20.csv")
)
# ## =========================================================
# ## 6. pseudobulk 생성
# ## =========================================================
# DefaultAssay(seu_pb) <- assay_use
# seu_pb$.pb_sample <- as.character(seu_pb@meta.data[[sample_col]])
# 
# pb_list <- AggregateExpression(
#   object = seu_pb,
#   assays = assay_use,
#   group.by = ".pb_sample",
#   return.seurat = FALSE,
#   verbose = FALSE
# )
# 
# pb_counts <- as.matrix(pb_list[[assay_use]])
# 
# cat("\nraw pseudobulk matrix dim (genes x samples):\n")
# print(dim(pb_counts))
# 
# cat("\nraw pseudobulk sample names:\n")
# print(colnames(pb_counts))
# keep_samples <- gsub('_','-',keep_samples)
# pb_counts <- pb_counts[, keep_samples[keep_samples %in% colnames(pb_counts)], drop = FALSE]
# 
# cat("\nreordered pseudobulk dim:\n")
# print(dim(pb_counts))
# 
# ## =========================================================
# ## 7. gene symbol uppercase + duplicate collapse
# ## =========================================================
# rownames(pb_counts) <- toupper(rownames(pb_counts))
# pb_counts2 <- rowsum(pb_counts, group = rownames(pb_counts), reorder = FALSE)
# pb_counts2 <- as.matrix(pb_counts2)
# 
# cat("\nuppercase/collapsed pseudobulk dim:\n")
# print(dim(pb_counts2))
# 
# saveRDS(pb_counts2, file.path(outdir, "04_pseudobulk_counts.rds"))
# write.csv(
#   pb_counts2[1:min(20, nrow(pb_counts2)), , drop = FALSE],
#   file.path(outdir, "04_pseudobulk_counts_head20.csv")
# )
# 
# ## =========================================================
# ## 8. logCPM 생성
# ## =========================================================
# pb_logcpm <- edgeR::cpm(pb_counts2, log = TRUE, prior.count = 1)
# 
# cat("\nlogCPM dim:\n")
# print(dim(pb_logcpm))
# 
# saveRDS(pb_logcpm, file.path(outdir, "05_pseudobulk_logcpm.rds"))
# write.csv(
#   pb_logcpm[1:min(20, nrow(pb_logcpm)), , drop = FALSE],
#   file.path(outdir, "05_pseudobulk_logcpm_head20.csv")
# )

## =========================================================
## 9. gene set 준비
##    - 5 hallmark
##    - 4 APC blocks
## =========================================================
hallmark_targets <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

hallmark_df <- msigdbr(
  species = species_use,
  collection = "H"
) %>%
  filter(gs_name %in% hallmark_targets) %>%
  transmute(
    gs_name,
    gene_symbol = toupper(gene_symbol)
  ) %>%
  distinct()

hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
hallmark_sets <- hallmark_sets[hallmark_targets]

ap_mhcii_core <- c(
  "CD74",
  "HLA-DMA", "HLA-DMB",
  "HLA-DOA", "HLA-DOB",
  "HLA-DPA1", "HLA-DPB1",
  "HLA-DQA1", "HLA-DQA2",
  "HLA-DQB1", "HLA-DQB2",
  "HLA-DRA",
  "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "CD74", "CIITA",
  "HLA-DMA", "HLA-DMB",
  "HLA-DOA", "HLA-DOB",
  "IFI30", "LGMN",
  "CTSB", "CTSD", "CTSL", "CTSS"
)

ap_ifn_immunoproteasome <- c(
  "TAP1", "TAP2", "B2M", "PSMB8", "PSMB9", "HLA-A", "HLA-B", "HLA-C"
)

ap_costim_activation <- c(
  "CD80", "CD86", "CD40", "ICAM1"
)

custom_sets <- list(
  AP_MHCII_CORE = ap_mhcii_core,
  #AP_MHCII_PROCESS = ap_mhcii_process,
  AP_IFN_IMMUNOPROTEASOME = ap_ifn_immunoproteasome,
  AP_COSTIM_ACTIVATION = ap_costim_activation
)

custom_sets <- lapply(custom_sets, toupper)

gene_sets <- c(hallmark_sets, custom_sets)

cat("\n전체 gene set 이름:\n")
print(names(gene_sets))
db_species_use <- 'HS'
go_collection <- 'C5'
reactome_collection <- 'C2'
bcell_go_targets <- c(
  "GOBP_B_CELL_ACTIVATION",
  "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_B_CELL_PROLIFERATION",
  "GOBP_GERMINAL_CENTER_FORMATION",
  "GOBP_SOMATIC_DIVERSIFICATION_OF_IMMUNOGLOBULINS",
  "GOBP_SOMATIC_DIVERSIFICATION_OF_IMMUNE_RECEPTORS_VIA_SOMATIC_MUTATION",
  "GOBP_PLASMA_CELL_DIFFERENTIATION",
  "GOBP_IMMUNOGLOBULIN_PRODUCTION",
  "GOBP_IMMUNOGLOBULIN_PRODUCTION_INVOLVED_IN_IMMUNOGLOBULIN_MEDIATED_IMMUNE_RESPONSE"
)

bcell_reactome_targets <- c(
  "REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR"
)

bcell_targets <- c(bcell_go_targets, bcell_reactome_targets)

norm_symbol <- function(x) {
  toupper(trimws(as.character(x)))
}
## GO BP gene sets
bcell_go_df <- msigdbr(
  db_species = db_species_use,
  species = species_use,
  collection = go_collection
) %>%
  filter(gs_name %in% bcell_go_targets) %>%
  transmute(
    gs_name,
    gene_symbol = norm_symbol(gene_symbol)
  ) %>%
  distinct()

## Reactome gene sets
bcell_reactome_df <- msigdbr(
  db_species = db_species_use,
  species = species_use,
  collection = reactome_collection
) %>%
  filter(gs_name %in% bcell_reactome_targets) %>%
  transmute(
    gs_name,
    gene_symbol = norm_symbol(gene_symbol)
  ) %>%
  distinct()

bcell_df <- bind_rows(bcell_go_df, bcell_reactome_df) %>%
  distinct()

## 실제로 불러온 set 확인
found_sets <- sort(unique(bcell_df$gs_name))
missing_sets <- setdiff(bcell_targets, found_sets)

cat("\n불러온 B-cell gene set:\n")
print(found_sets)

if (length(missing_sets) > 0) {
  cat("\n주의: msigdbr에서 못 찾은 gene set:\n")
  print(missing_sets)
}

gene_sets2 <- split(bcell_df$gene_symbol, bcell_df$gs_name)
gene_sets2 <- gene_sets2[intersect(bcell_targets, names(gene_sets2))]

cat("\n전체 gene set 이름:\n")
print(names(gene_sets))
gene_sets3 <- c(gene_sets, gene_sets2)
## =========================================================
## 10. gene set overlap 확인
## =========================================================
expr_genes <- rownames(pb_logcpm)

overlap_tbl <- tibble(
  gs_name = names(gene_sets3),
  n_input_genes = vapply(gene_sets3, length, integer(1)),
  n_overlap = vapply(gene_sets3, function(gs) sum(gs %in% expr_genes), integer(1))
) %>%
  arrange(desc(n_overlap))

cat("\ngene set overlap:\n")
print(overlap_tbl)
# gs_name                                                                            n_input_genes n_overlap
# <chr>                                                                                      <int>     <int>
#   1 GOBP_B_CELL_ACTIVATION                                                                       297       218
# 2 HALLMARK_INTERFERON_GAMMA_RESPONSE                                                           200       169
# 3 HALLMARK_TNFA_SIGNALING_VIA_NFKB                                                             200       139
# 4 REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR                                                153       136
# 5 HALLMARK_INFLAMMATORY_RESPONSE                                                               200       100
# 6 GOBP_IMMUNOGLOBULIN_PRODUCTION                                                               112        87
# 7 HALLMARK_INTERFERON_ALPHA_RESPONSE                                                            97        86
# 8 GOBP_B_CELL_PROLIFERATION                                                                    100        75
# 9 GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY                                                        79        68
# 10 GOBP_SOMATIC_DIVERSIFICATION_OF_IMMUNOGLOBULINS                                               70        58
# 11 HALLMARK_IL6_JAK_STAT3_SIGNALING                                                              87        54
# 12 GOBP_IMMUNOGLOBULIN_PRODUCTION_INVOLVED_IN_IMMUNOGLOBULIN_MEDIATED_IMMUNE_RESPONSE            59        46
# 13 AP_MHCII_CORE                                                                                 28        25
# 14 GOBP_SOMATIC_DIVERSIFICATION_OF_IMMUNE_RECEPTORS_VIA_SOMATIC_MUTATION                         16        13
# 15 GOBP_GERMINAL_CENTER_FORMATION                                                                14        12
# 16 AP_IFN_IMMUNOPROTEASOME                                                                        8         8
# 17 GOBP_PLASMA_CELL_DIFFERENTIATION                                                              10         7
# 18 AP_COSTIM_ACTIVATION                                                                           4         4
write.csv(overlap_tbl, file.path(outdir, "06_gene_set_overlap.csv"), row.names = FALSE)

keep_gs <- overlap_tbl %>%
  filter(n_overlap >= min_geneset_overlap) %>%
  pull(gs_name)

gene_sets_filt <- gene_sets3[keep_gs] %>%
  lapply(function(gs) intersect(gs, expr_genes))

saveRDS(gene_sets_filt, file.path(outdir, "06_gene_sets_filtered.rds"))

cat("\n최종 사용 gene set:\n")
print(names(gene_sets_filt))

## =========================================================
## 11. GMT 파일 저장
## =========================================================
gmt_file <- file.path(outdir, "blineage_sets.gmt")

gmt_lines <- vapply(names(gene_sets_filt), function(gs) {
  paste(c(gs, "na", gene_sets_filt[[gs]]), collapse = "\t")
}, character(1))

writeLines(gmt_lines, con = gmt_file)

cat("\nGMT file exists:\n")
print(file.exists(gmt_file))
print(gmt_file)

## =========================================================
## 12. GCT 파일 저장
## =========================================================
gct_file <- file.path(outdir, "blineage_expr.gct")

con <- file(gct_file, open = "wt")
writeLines("#1.2", con)
writeLines(sprintf("%d\t%d", nrow(pb_logcpm), ncol(pb_logcpm)), con)
writeLines(
  paste(c("Name", "Description", colnames(pb_logcpm)), collapse = "\t"),
  con
)

gct_df <- data.frame(
  Name = rownames(pb_logcpm),
  Description = "na",
  pb_logcpm,
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

cat("\nGCT file exists:\n")
print(file.exists(gct_file))
print(gct_file)

cat("\nGCT header preview:\n")
print(readLines(gct_file, n = 5))

## =========================================================
## 13. ssGSEA2 실행
## =========================================================
ssgsea_out_prefix <- "blineage_hallmark_ap"
ssgsea_log_file <- file.path(outdir, paste0(ssgsea_out_prefix, ".run.log"))

res_ssgsea <- ssGSEA2::run_ssGSEA2(
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

## sample 이름
sample_names <- colnames(pb_logcpm)

## pathway별 nested list에 sample names 부여
res_ssgsea_named <- lapply(res_ssgsea, function(pathway_obj) {
  n_use <- min(length(pathway_obj), length(sample_names))
  names(pathway_obj)[seq_len(n_use)] <- sample_names[seq_len(n_use)]
  pathway_obj
})

## 확인
names(res_ssgsea_named[["AP_MHCII_CORE"]])[1:5]

saveRDS(res_ssgsea, file.path(outdir, "07_ssgsea2_raw_result.rds"))

cat("\nssGSEA2 result structure:\n")
str(res_ssgsea, max.level = 2)

cat("\nssGSEA2 log file exists:\n")
print(file.exists(ssgsea_log_file))

ssgsea_long <- bind_rows(lapply(names(res_ssgsea_named), function(gs_name) {
  pathway_obj <- res_ssgsea_named[[gs_name]]
  
  tibble(
    sample = names(pathway_obj),
    gs_name = gs_name,
    ssgsea2_es = vapply(pathway_obj, function(x) {
      if (!is.null(x$ES)) as.numeric(x$ES[1]) else NA_real_
    }, numeric(1))
  )
}))

head(ssgsea_long)

## =========================================================
## 14. singscore 계산
## =========================================================
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
  
  tibble(
    sample = rownames(ss),
    gs_name = gs_name,
    singscore = ss[[score_col]]
  )
})

singscore_long <- bind_rows(singscore_long_list)

cat("\nsingscore long preview:\n")
print(head(singscore_long))

write.csv(singscore_long, file.path(outdir, "08_singscore_long.csv"), row.names = FALSE)

## =========================================================
## 15. PC1 / PC2 / eigengene 계산
## =========================================================
pc_long_list <- lapply(names(gene_sets_filt), function(gs_name) {
  genes <- intersect(gene_sets_filt[[gs_name]], rownames(pb_logcpm))
  X <- t(pb_logcpm[genes, , drop = FALSE])   # sample x genes
  
  if (ncol(X) == 0) {
    return(tibble(
      sample = rownames(X),
      gs_name = gs_name,
      pc1 = NA_real_,
      pc2 = NA_real_,
      eigengene = NA_real_
    ))
  }
  
  if (nrow(X) < 2) {
    pc1 <- rep(0, nrow(X))
    pc2 <- rep(NA_real_, nrow(X))
    names(pc1) <- rownames(X)
  } else if (ncol(X) == 1) {
    pc1 <- as.numeric(scale(X[, 1]))
    pc2 <- rep(NA_real_, length(pc1))
    names(pc1) <- rownames(X)
  } else {
    pca <- prcomp(X, center = TRUE, scale. = TRUE)
    pc1 <- pca$x[, 1]
    pc2 <- if (ncol(pca$x) >= 2) pca$x[, 2] else rep(NA_real_, nrow(X))
    
    avg_expr <- rowMeans(X, na.rm = TRUE)
    cc <- suppressWarnings(cor(pc1, avg_expr, use = "pairwise.complete.obs"))
    if (!is.na(cc) && cc < 0) {
      pc1 <- -pc1
    }
  }
  
  tibble(
    sample = rownames(X),
    gs_name = gs_name,
    pc1 = as.numeric(pc1),
    pc2 = as.numeric(pc2),
    eigengene = as.numeric(pc1)
  )
})

pc_long <- bind_rows(pc_long_list)

cat("\nPC score long preview:\n")
print(head(pc_long))

write.csv(pc_long, file.path(outdir, "09_pc_long.csv"), row.names = FALSE)

## =========================================================
## 16. partial merge
## =========================================================
feature_long <- ssgsea_long %>%
  full_join(singscore_long, by = c("sample", "gs_name")) %>%
  full_join(pc_long, by = c("sample", "gs_name")) %>%
  left_join(overlap_tbl, by = "gs_name")

write.csv(feature_partial, file.path(outdir, "10_feature_partial_no_ssgsea.csv"), row.names = FALSE)
write.csv(feature_long, file.path(outdir, "10_feature_partial_yes_ssgsea.csv"), row.names = FALSE)

cat("\n중간 산물 저장 완료:\n")
print(list.files(outdir))
########################################
# 2. de novo gene set
# NMF
##################################3
## =========================================================
## 0. packages
## =========================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(GeneNMF)
library(edgeR)
library(singscore)
library(RColorBrewer)
library(patchwork)

## =========================================================
## 1. input
## =========================================================
seu_use <- b_subust_clean   # clean B-lineage Seurat object
sample_col <- "sample"
cohort_col <- "cohort_l2"
target_cohort <- "Lee_p1_base"
assay_use <- "RNA"

outdir <- "data/20260309_pilot/results/version2/B/geneNMF_stepwise"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_sample <- 20

## GeneNMF 설정
k_vec <- 4:7
nfeatures_hvg <- 2000
nmf_seed <- 123

## B-cell에서 제외 권장 gene blocklist
## - technical / clone-specific bias 완화 목적
## - 필요시 더 추가
gene_blocklist_patterns <- c(
  "^MT-",
  "^RPS", "^RPL",
  "^HSP",
  "^IG[HKL]",
  "^TRA", "^TRB", "^TRD", "^TRG"
)

## =========================================================
## 2. cohort filtering
## =========================================================
meta0 <- seu_use@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

cat("cohort distribution:\n")
print(table(meta0[[cohort_col]], useNA = "ifany"))

keep_cells_cohort <- meta0$cell[
  !is.na(meta0[[cohort_col]]) & meta0[[cohort_col]] == target_cohort
]

seu_cohort <- subset(seu_use, cells = keep_cells_cohort)

cat("\nCells in target cohort:\n")
print(ncol(seu_cohort))

## =========================================================
## 3. sample filtering
## =========================================================
meta1 <- seu_cohort@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

cat("\nsample distribution in target cohort:\n")
print(table(meta1[[sample_col]], useNA = "ifany"))

keep_cells_non_na <- meta1$cell[!is.na(meta1[[sample_col]]) & meta1[[sample_col]] != ""]
seu_non_na <- subset(seu_cohort, cells = keep_cells_non_na)

meta2 <- seu_non_na@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

sample_meta2 <- meta2 %>%
  transmute(sample = as.character(.data[[sample_col]])) %>%
  count(sample, name = "n_cells") %>%
  arrange(desc(n_cells))

cat("\nSample counts after NA removal:\n")
print(sample_meta2)

keep_samples <- sample_meta2 %>%
  filter(n_cells >= min_cells_per_sample) %>%
  pull(sample)

cat("\nSamples passing min_cells_per_sample:\n")
print(keep_samples)

if (length(keep_samples) < 2) {
  stop("Need at least 2 samples for robust multi-sample GeneNMF consensus analysis.")
}

keep_cells <- meta2$cell[meta2[[sample_col]] %in% keep_samples]
seu_nm <- subset(seu_non_na, cells = keep_cells)

## =========================================================
## 4. split into per-sample Seurat objects
## =========================================================
obj_list <- SplitObject(seu_nm, split.by = sample_col)

## sample 이름 확인
cat("\nObject list names:\n")
print(names(obj_list))

## 혹시 sample별 cell 수 재확인
sample_sizes <- sapply(obj_list, ncol)
cat("\nCell counts per sample object:\n")
print(sample_sizes)

## =========================================================
## 5. optional: keep only reasonably sized samples
##    (already filtered, but one more explicit check)
## =========================================================
obj_list <- obj_list[sample_sizes >= min_cells_per_sample]

cat("\nFinal sample objects for GeneNMF:\n")
print(names(obj_list))

## =========================================================
## 6. Gene blocklist 만들기
## =========================================================
all_genes <- rownames(obj_list[[1]])
genes_blocklist <- unique(unlist(lapply(gene_blocklist_patterns, function(p) {
  grep(p, all_genes, value = TRUE, ignore.case = FALSE)
})))

cat("\nNumber of blocked genes:\n")
print(length(genes_blocklist))
write.csv(
  data.frame(gene = genes_blocklist),
  file.path(outdir, "01_genes_blocklist.csv"),
  row.names = FALSE
)

## =========================================================
## 7a. Seurat v5 multi-layer assay 정리
##     GeneNMF가 GetAssayData()로 읽을 수 있게 RNA assay를 single-layer로 만듦
## =========================================================
obj_list <- lapply(obj_list, function(x) {
  DefaultAssay(x) <- assay_use
  
  cat("\nBefore JoinLayers:\n")
  print(Layers(x[[assay_use]]))
  
  # 여러 layer가 있으면 합치기
  x <- JoinLayers(x, assay = assay_use)
  
  cat("After JoinLayers:\n")
  print(Layers(x[[assay_use]]))
  
  # data layer가 없으면 다시 만들어줌
  if (!"data" %in% Layers(x[[assay_use]])) {
    x <- NormalizeData(x, assay = assay_use, verbose = FALSE)
  }
  
  x
})
## =========================================================
## 7. HVG filtering with GeneNMF helper
##    each sample별로 variable genes 계산
## =========================================================

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

## sample별 HVG 확인
hvg_list <- lapply(obj_list_hvg, VariableFeatures)
hvg_summary <- data.frame(
  sample = names(hvg_list),
  n_hvg = sapply(hvg_list, length)
)

cat("\nHVG summary:\n")
print(hvg_summary)
write.csv(hvg_summary, file.path(outdir, "02_hvg_summary.csv"), row.names = FALSE)

## 공통/합집합 HVG도 확인
hvg_union <- sort(unique(unlist(hvg_list)))
cat("\nUnion HVG size:\n")
print(length(hvg_union))

write.csv(
  data.frame(gene = hvg_union),
  file.path(outdir, "02_hvg_union.csv"),
  row.names = FALSE
)

## =========================================================
## 8. run multiNMF
## =========================================================
## GeneNMF docs:
## multiNMF(obj.list, assay="RNA", slot="data", k=5:6, hvg=NULL, nfeatures=2000,
##          L1=c(0,0), min.exp=0.01, max.exp=3, center=FALSE, scale=FALSE,
##          min.cells.per.sample=10, hvg.blocklist=NULL, seed=123)
## multiNMF는 sample별로 여러 k에서 NMF를 돌립니다.
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

cat("\nGeneNMF programs object saved.\n")
cat("Top-level names:\n")
print(names(geneNMF_programs))

## =========================================================
## 9. extract top genes per NMF program
## =========================================================
## getNMFgenes(): top genes for each individual program
nmf_genes <- GeneNMF::getNMFgenes(
  nmf.res = geneNMF_programs,
  specificity.weight = 5,
  weight.explained = 0.5,
  max.genes = 100
)

saveRDS(nmf_genes, file.path(outdir, "04_nmf_genes.rds"))

cat("\nNumber of individual programs:\n")
print(length(nmf_genes))

## preview
print(utils::head(nmf_genes, 3))

## program gene list를 csv로 풀어 저장
nmf_genes_long <- bind_rows(lapply(names(nmf_genes), function(prog) {
  data.frame(
    program = prog,
    gene = nmf_genes[[prog]],
    stringsAsFactors = FALSE
  )
}))
write.csv(nmf_genes_long, file.path(outdir, "04_nmf_genes_long.csv"), row.names = FALSE)

## =========================================================
## 10. derive consensus metaprograms
## =========================================================
## getMetaPrograms()는 여러 샘플 / 여러 k에 걸쳐 recurrent program을 묶습니다.
geneNMF_metaprograms <- GeneNMF::getMetaPrograms(
  nmf.res = geneNMF_programs,
  nMP = 6,                  # 일단 시작값; 필요시 4~10 비교
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

anno_colors <- brewer.pal(n=6, name="Paired")
names(anno_colors) <- names(geneNMF_metaprograms$metaprograms.genes)
ph <- plotMetaPrograms(geneNMF_metaprograms, annotation_colors = anno_colors)
ph #8x9 #20260325_B_NMF_metaprogam_heatmap
saveRDS(geneNMF_metaprograms, file.path(outdir, "05_geneNMF_metaprograms.rds"))

cat("\nMeta-program object names:\n")
print(names(geneNMF_metaprograms))

## meta-program genes
mp_genes <- geneNMF_metaprograms$metaprograms.genes
mp_metrics <- geneNMF_metaprograms$metaprograms.metrics

cat("\nMeta-program names:\n")
print(names(mp_genes))

cat("\nMeta-program metrics:\n")
print(mp_metrics)

write.csv(mp_metrics, file.path(outdir, "05_metaprogram_metrics.csv"), row.names = FALSE)

mp_genes_long <- bind_rows(lapply(names(mp_genes), function(mp) {
  data.frame(
    metaprogram = mp,
    gene = mp_genes[[mp]],
    stringsAsFactors = FALSE
  )
}))
write.csv(mp_genes_long, file.path(outdir, "05_metaprogram_genes_long.csv"), row.names = FALSE)

## =========================================================
## 11. optional: visualize meta-program similarity
## =========================================================
## plot object 저장
# pdf(file.path(outdir, "06_metaprogram_heatmap.pdf"), width = 8, height = 7)
# print(GeneNMF::plotMetaPrograms(geneNMF_metaprograms))
# dev.off()

## =========================================================
## 12. build pseudobulk matrix for the same sample set
##    (for patient-level feature scoring)
## =========================================================
DefaultAssay(seu_nm) <- assay_use

pb_list <- AggregateExpression(
  object = seu_nm,
  assays = assay_use,
  group.by = sample_col,
  return.seurat = FALSE,
  verbose = FALSE
)

pb_counts <- as.matrix(pb_list[[assay_use]])
keep_samples <- gsub('_','-',keep_samples)
pb_counts <- pb_counts[, keep_samples[keep_samples %in% colnames(pb_counts)], drop = FALSE]

rownames(pb_counts) <- toupper(rownames(pb_counts))
pb_counts <- rowsum(pb_counts, group = rownames(pb_counts), reorder = FALSE)
pb_counts <- as.matrix(pb_counts)

pb_logcpm <- edgeR::cpm(pb_counts, log = TRUE, prior.count = 1)

saveRDS(pb_counts, file.path(outdir, "07_pseudobulk_counts.rds"))
saveRDS(pb_logcpm, file.path(outdir, "07_pseudobulk_logcpm.rds"))

cat("\nPseudobulk matrix dim:\n")
print(dim(pb_logcpm))

## =========================================================
## 13. score metaprograms on pseudobulk using singscore
##    (sample-level features)
## =========================================================
rank_data <- singscore::rankGenes(pb_logcpm, tiesMethod = "min")

mp_singscore_long <- bind_rows(lapply(names(mp_genes), function(mp) {
  genes <- intersect(toupper(mp_genes[[mp]]), rownames(pb_logcpm))
  if (length(genes) < 3) return(NULL)
  
  ss <- singscore::simpleScore(
    rankData = rank_data,
    upSet = genes,
    knownDirection = TRUE
  )
  
  score_col <- if ("TotalScore" %in% colnames(ss)) {
    "TotalScore"
  } else if ("UpScore" %in% colnames(ss)) {
    "UpScore"
  } else {
    colnames(ss)[1]
  }
  
  tibble(
    sample = rownames(ss),
    metaprogram = mp,
    singscore = ss[[score_col]],
    n_overlap = length(genes)
  )
}))

write.csv(mp_singscore_long, file.path(outdir, "08_metaprogram_singscore_long.csv"), row.names = FALSE)

## wide feature table
mp_feature_wide <- mp_singscore_long %>%
  mutate(mp_key = tolower(gsub("[^A-Za-z0-9]+", "_", metaprogram))) %>%
  select(sample, mp_key, singscore, n_overlap) %>%
  pivot_wider(
    names_from = mp_key,
    values_from = c(singscore, n_overlap),
    names_glue = "{mp_key}__{.value}"
  )

write.csv(mp_feature_wide, file.path(outdir, "08_metaprogram_feature_wide.csv"), row.names = FALSE)

cat("\nFinal GeneNMF-based sample feature table preview:\n")
print(head(mp_feature_wide))

####################################
# MERGE FEATURES
####################################
# res: calc_b_continuum_score
# feature_partial: singscore, epigene based score, (ssGSEA는 없)
# NMF: mp_feature_wide/mp_singscore_long

library(dplyr)
library(tidyr)
library(tibble)

## -------------------------------------------------------
## 0. helper
## -------------------------------------------------------
sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

## -------------------------------------------------------
## 1. continuum score table
##    calc_b_continuum_score_v2 결과에서 sample당 1행 추출
##    >> continuum_df
## -------------------------------------------------------
## 예: res <- calc_b_continuum_score_v2(...)
continuum_df <- res$patient_summary %>%
  rename(sample = .sample)

cat("continuum_df dim:\n")
print(dim(continuum_df))
print(head(continuum_df))
colnames(continuum_df) <- paste0('Centroid_',colnames(continuum_df))
colnames(continuum_df)[1] <- 'sample'
## -------------------------------------------------------
## 2. curated gene set scores
##    A) 이미 wide면 그대로 사용
##    B) long이면 wide로 변환
##    >> curated_feature_wide_full
## -------------------------------------------------------

## ---- 경우 A: score_wide가 이미 있으면 ----
## curated_feature_wide <- res_prog$score_wide

## ---- 경우 B1: feature_partial만 있고, 여기엔 singscore/pc1/pc2/eigengene만 있을 때 ----
## expected columns: sample, gs_name, singscore, pc1, pc2, eigengene
if (exists("feature_partial")) {
  curated_feature_wide_partial <- feature_partial %>%
    mutate(gs_key = sanitize_name(gs_name)) %>%
    select(sample, gs_key, singscore, pc1, pc2, eigengene) %>%
    pivot_wider(
      names_from = gs_key,
      values_from = c(singscore, pc1, pc2, eigengene),
      names_glue = "{gs_key}__{.value}"
    )
} else {
  curated_feature_wide_partial <- NULL
}
## ---- 경우 B2: curated long table에 ssGSEA/singscore/pc1/pc2/eigengene이 다 있을 때 ----
## expected columns:
## sample, gs_name, ssgsea2_nes, singscore, pc1, pc2, eigengene
if (exists("feature_long")) {
  curated_feature_wide_full <- feature_long %>%
    mutate(gs_key = sanitize_name(gs_name)) %>%
    select(sample, gs_key, ssgsea2_es, singscore, pc1, pc2, eigengene) %>%
    pivot_wider(
      names_from = gs_key,
      values_from = c(ssgsea2_es, singscore, pc1, pc2, eigengene),
      names_glue = "{gs_key}__{.value}"
    )
} else {
  curated_feature_wide_full <- NULL
}
colnames(curated_feature_wide_full) <- paste0('curated_gene_',colnames(curated_feature_wide_full))
colnames(curated_feature_wide_full)[1] <- 'sample'

## 우선순위: full > partial
curated_feature_wide <- if (!is.null(curated_feature_wide_full)) {
  curated_feature_wide_full
} else if (!is.null(curated_feature_wide_partial)) {
  curated_feature_wide_partial
} else {
  NULL
}

if (!is.null(curated_feature_wide)) {
  cat("curated_feature_wide dim:\n")
  print(dim(curated_feature_wide))
  print(head(curated_feature_wide[, 1:min(8, ncol(curated_feature_wide)), drop = FALSE]))
}
## -------------------------------------------------------
## 3. GeneNMF module scores
##    A) 이미 wide면 그대로 사용
##    B) long이면 wide로 변환
##    >> mp_feature_wide
## -------------------------------------------------------

## ---- 경우 A: mp_feature_wide가 이미 있으면 ----
## nmf_feature_wide <- mp_feature_wide
colnames(mp_feature_wide) <- paste0('denovo_gene_',colnames(mp_feature_wide))
colnames(mp_feature_wide)[1] <-'sample'
## ---- 경우 B: mp_singscore_long 같은 long table만 있으면 ----
## expected columns: sample, metaprogram, singscore, n_overlap
if (exists("mp_feature_wide")) {
  nmf_feature_wide <- mp_feature_wide
} else if (exists("mp_singscore_long")) {
  nmf_feature_wide <- mp_singscore_long %>%
    mutate(mp_key = sanitize_name(metaprogram)) %>%
    select(sample, mp_key, singscore, n_overlap) %>%
    pivot_wider(
      names_from = mp_key,
      values_from = c(singscore, n_overlap),
      names_glue = "nmf_{mp_key}__{.value}"
    )
} else {
  nmf_feature_wide <- NULL
}

if (!is.null(nmf_feature_wide)) {
  cat("nmf_feature_wide dim:\n")
  print(dim(nmf_feature_wide))
  print(head(nmf_feature_wide[, 1:min(8, ncol(nmf_feature_wide)), drop = FALSE]))
}

## -------------------------------------------------------
## 4. sample consistency check
## -------------------------------------------------------
continuum_df$sample <- gsub('_','-',continuum_df$sample)
sample_sets <- list(
  continuum = continuum_df$sample,
  curated = if (!is.null(curated_feature_wide_full)) curated_feature_wide_full$sample else character(0),
  nmf = if (!is.null(mp_feature_wide)) mp_feature_wide$sample else character(0)
)

cat("\nSample overlap check:\n")
print(lapply(sample_sets, unique))

if (!is.null(curated_feature_wide_full)) {
  cat("\nSamples in continuum but not curated:\n")
  print(setdiff(continuum_df$sample, curated_feature_wide_full$sample))
  cat("\nSamples in curated but not continuum:\n")
  print(setdiff(curated_feature_wide_full$sample, continuum_df$sample))
}

if (!is.null(mp_feature_wide)) {
  cat("\nSamples in continuum but not NMF:\n")
  print(setdiff(continuum_df$sample, mp_feature_wide$sample))
  cat("\nSamples in NMF but not continuum:\n")
  print(setdiff(mp_feature_wide$sample, continuum_df$sample))
}
## -------------------------------------------------------
## 5. merge
## -------------------------------------------------------
final_feature_df <- continuum_df

if (!is.null(curated_feature_wide_full)) {
  final_feature_df <- final_feature_df %>%
    left_join(curated_feature_wide_full, by = "sample")
}

if (!is.null(mp_feature_wide)) {
  final_feature_df <- final_feature_df %>%
    left_join(mp_feature_wide, by = "sample")
}

cat("\nfinal_feature_df dim:\n")
print(dim(final_feature_df)) #79x54
print(head(final_feature_df[, 1:min(12, ncol(final_feature_df)), drop = FALSE]))

write.csv(
  final_feature_df,
  "data/20260309_pilot/results/version2/B/final_feature_df_merged.csv",
  row.names = FALSE
)
 
#################
# data for CellRank
#################

ggplot(b_sbust@meta.data, aes(x=sample,fill=manual.cluster_l2))+geom_bar(position='fill')+theme(axis.text.x= element_text(angle=90,hjust=1,vjust=1,size=3))#4x8 #20260325_B_subset_composition
library(Matrix)
OUTDIR <- "data/20260309_pilot/results/version2/B/CellRank_input/"
ASSAY_USE <- "RNA"

# metadata column names: 네 object에 맞게 수정
COHORT_COL <- "cohort_l2"
TARGET_COHORT <- "Lee_p1_base"

SAMPLE_COL <- "sample"              # patient ID column
STATE_COL  <- "manual.cluster_l2"       # naive / memory / intermediate_Resting

KEEP_STATES <- c("naive", "memory", "intermediate_resting","plasma")
MIN_CELLS_PER_PATIENT <- 20

# =========================================================
# 1. cohort + clean B states subset
# =========================================================
md <- b_sbust@meta.data

keep_cells <- rownames(md)[
  md[[COHORT_COL]] == TARGET_COHORT &
    md[[STATE_COL]] %in% KEEP_STATES
]

b_sbust_b <- subset(b_sbust, cells = keep_cells)

# =========================================================
# 2. patient filter: at least 20 cells in the retained subset
# =========================================================
patient_n <- as.data.table(b_sbust_b@meta.data)[, .N, by = SAMPLE_COL]
setnames(patient_n, "N", "n_cells_retained")

eligible_patients <- patient_n[n_cells_retained >= MIN_CELLS_PER_PATIENT][[SAMPLE_COL]]

b_sbust_b <- subset(
  b_sbust_b,
  cells = rownames(b_sbust_b@meta.data)[b_sbust_b@meta.data[[SAMPLE_COL]] %in% eligible_patients]
)

# re-check
patient_n2 <- as.data.table(b_sbust_b@meta.data)[, .N, by = SAMPLE_COL]
setnames(patient_n2, "N", "n_cells_retained")
fwrite(patient_n2, file.path(OUTDIR, "patient_counts_after_filter.csv"))

# =========================================================
# 3. export counts + obs + var
#    Python에서 AnnData로 다시 조립할 예정
# =========================================================
counts <- GetAssayData(b_sbust_b, assay = ASSAY_USE, layer = "counts")

# 너무 희귀한 gene은 제거
gene_keep <- Matrix::rowSums(counts > 0) >= 3
counts <- counts[gene_keep, ]

obs <- b_sbust_b@meta.data
obs$cell_id <- rownames(obs)
obs <- obs[colnames(counts), , drop = FALSE]

# var table
var <- data.table(
  gene_id = rownames(counts),
  gene_symbol = rownames(counts)
)

# write
writeMM(counts, file.path(OUTDIR, "counts.mtx"))
fwrite(as.data.table(obs), file.path(OUTDIR, "obs.csv"))
fwrite(var, file.path(OUTDIR, "var.csv"))

# optional: keep cell IDs in separate file too
fwrite(data.table(cell_id = colnames(counts)), file.path(OUTDIR, "barcodes.csv"))

# save tiny metadata summary
summary_dt <- data.table(
  n_cells = ncol(counts),
  n_genes = nrow(counts),
  n_patients = length(unique(obs[[SAMPLE_COL]])),
  cohort = TARGET_COHORT
)
fwrite(summary_dt, file.path(OUTDIR, "export_summary.csv"))

cat("Done: exported clean B subset for CellRank\n")

##############################
## Merge CellRank to final_feature_df
##############################
# final_feature_df
# CellRank patient-level output
CELLRANK_PATIENT_CSV <- "data/20260309_pilot/results/version2/B/CellRank_stepwise/cellrank_b3_patient_features.csv"

OUT_RDS <- "data/20260309_pilot/results/version2/B/CellRank_stepwise/final_feature_df_with_cellrank_b3.rds"
OUT_CSV <- "data/20260309_pilot/results/version2/B/CellRank_stepwise/final_feature_df_with_cellrank_b3.csv"

cellrank_df <- fread(CELLRANK_PATIENT_CSV)
cellrank_df$sample <- gsub('_','-',cellrank_df$sample)
cat("Loaded final_feature_df:", nrow(final_feature_df), "rows x", ncol(final_feature_df), "cols\n")
cat("Loaded cellrank_df    :", nrow(cellrank_df), "rows x", ncol(cellrank_df), "cols\n")

patient_col='sample'
# CellRank 컬럼만 추리기
cellrank_df <- cellrank_df %>%
  select(
    all_of(patient_col),
    starts_with("cr_b3_"),
    starts_with("qc_b3_")
  ) %>%
  distinct(.data[[patient_col]], .keep_all = TRUE)

colnames(cellrank_df) <- paste0('cellrank_',colnames(cellrank_df))
colnames(cellrank_df)[1] <- 'sample'


# 충돌 컬럼 체크
overlap_cols <- intersect(
  setdiff(colnames(final_feature_df), patient_col),
  setdiff(colnames(cellrank_df), patient_col)
)


final_feature_df_cellrank <- final_feature_df %>%
  left_join(cellrank_df, by = patient_col)

saveRDS(final_feature_df_cellrank, OUT_RDS)
fwrite(as.data.table(final_feature_df_cellrank), OUT_CSV)

final_feature_df_cellrank_filt <- final_feature_df_cellrank[,!(colnames(final_feature_df_cellrank) %in% c('curated_gene_ap_ifn_immunoproteasome__pc1', 'curated_gene_ap_costim_activation__pc1', 'curated_gene_ap_ifn_immunoproteasome__pc2', 'curated_gene_ap_costim_activation__pc2', 'curated_gene_ap_ifn_immunoproteasome__eigengene', 'curated_gene_ap_costim_activation__eigengene', 'denovo_gene_mp1__n_overlap', 'denovo_gene_mp2__n_overlap', 'denovo_gene_mp3__n_overlap', 'Centroid_n_b_cells', 'Centroid_frac_intermediate_like', 'Centroid_observed_frac_intermediate', 'Centroid_keep_for_analysis', 'cellrank_cr_b3_priming_mean', 'cellrank_cr_b3_priming_q75', 'cellrank_cr_b3_n_cells', 'cellrank_qc_b3_subset_frac_intermediate_resting', 'cellrank_qc_b3_subset_frac_memory', 'cellrank_qc_b3_subset_frac_naive', 'cellrank_qc_b3_subset_frac_plasma'))]
OUT_CSV <- "data/20260309_pilot/results/version2/B/CellRank_stepwise/final_feature_df_with_cellrank_b3_filt.csv"
fwrite(as.data.table(final_feature_df_cellrank_filt), OUT_CSV)
