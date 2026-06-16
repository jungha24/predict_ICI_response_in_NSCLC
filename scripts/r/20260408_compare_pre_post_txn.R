BiocManager::install("miloR")
library(miloR)
library(Seurat)
library(SingleCellExperiment)
BiocManager::install("scater")
library(scater)
BiocManager::install("scran")
library(scran)
library(dplyr)
library(data.table)
#####
load('data/20260309_pilot/results/20260310_integrated.RData')
DimPlot(seu_integrated,reduction='umap',label=T)

meta <- fread('data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt')
meta$sample_id <- gsub('_','-',meta$sample_id)

# subset 파일들
subset_files <- c(
  "data/20260309_pilot/results/version2/B/b_sbust.rds",
  "data/20260309_pilot/results/version2/NK/NK_subset.rds",
  "data/20260309_pilot/results/version2/Monocyte/Mono_subset.rds",
  "data/20260309_pilot/results/version2/CD4T/CD4T_subset.rds",
  "data/20260309_pilot/results/version2/CD8T/CD8T_subset.rds",
  "data/20260309_pilot/results/version2/nonconventional_T/nonconventional_T_subset.rds"
)

anno_col <- "manual.cluster_l2"

anno_list <- lapply(subset_files, function(f) {
  obj <- readRDS(f)
  
  # subset object 안에서 cell 이름과 annotation 추출
  dt <- data.table(
    cell = Cells(obj),
    annotation = obj[[anno_col, drop = TRUE]]
  )
  dt
})

anno_dt <- rbindlist(anno_list, use.names = TRUE, fill = TRUE)

# 중복 cell 체크
dup_cells <- anno_dt[duplicated(cell), unique(cell)]
length(dup_cells)

anno_dt_good <- anno_dt[!(anno_dt$cell %in% dup_cells),]
anno_dt_dup <- anno_list[[6]]
anno_dt <- rbind(anno_dt_good, anno_dt_dup)


anno_df <- as.data.frame(anno_dt)
rownames(anno_df) <- anno_df$cell
anno_df$cell <- NULL
colnames(anno_df) <- anno_col

anno_col <- "manual.cluster_l2" 
new_colname <- "manual.cluster_l2"

# 1) seu_integrated의 기존 manual.cluster를 기본값으로 사용
full_vec <- setNames(
  as.character(seu_integrated$manual.cluster),
  Cells(seu_integrated)
)

# 2) anno_df에 있는 cell만 새 annotation으로 덮어쓰기
common_cells <- intersect(rownames(anno_df), Cells(seu_integrated))
full_vec[common_cells] <- as.character(anno_df[common_cells, anno_col])

# 3) metadata로 추가
seu_integrated <- AddMetaData(
  seu_integrated,
  metadata = full_vec,
  col.name = new_colname
)
DimPlot(seu_integrated_subset,reduction='umap',group.by='manual.cluster_l2',label=T)

seu_integrated_subset <- subset(seu_integrated, subset = !(manual.cluster_l2 %in% c('ambiguous','artifact','contam','platelet-like contam','Tcell contam','platelet')))#6x10 #20260410_manual_cluster_l2_wo_contam
seu_integrated_subset <- subset(seu_integrated_subset, subset = cohort_l2 != "AIDA_")
seu_integrated_subset <- subset(seu_integrated_subset, subset = !(sample %in% c("IO_SC_056","IO_SC_116","IO_SC_118","IO_SC_126","IO_SC_127","Normal")))

meta$sample_id <- gsub('-','_',meta$sample_id)
dcb_sample <- meta[meta$`Binarized response`=='DCB',]$sample_id
ncb_sample <- meta[meta$`Binarized response`=='NCB',]$sample_id
seu_integrated_subset$response <- ifelse(seu_integrated_subset$sample %in% dcb_sample, 'DCB',
                                         ifelse(seu_integrated_subset$sample %in% ncb_sample,'NCB','unknown'))


seu_integrated_subset <- subset(seu_integrated_subset, subset = seu_integrated_subset$response != 'unknown')

### miloR
sce <- as.SingleCellExperiment(seu_integrated_subset)

reducedDim(sce, "PCA") <- Embeddings(seu_integrated_subset, reduction = "pca")
reducedDim(sce, "UMAP") <- Embeddings(seu_integrated_subset, reduction = "umap")

milo <- Milo(sce)

milo <- buildGraph(
  milo,
  k = 20,
  d = 20,
  reduced.dim = "PCA"
)

milo <- makeNhoods(
  milo,
  prop = 0.1,
  k = 20,
  d = 20,
  reduced_dims = "PCA",
  refined = TRUE
)
saveRDS(milo, file='data/20260309_pilot/results/milo_20260408.rds')
milo
plotNhoodSizeHist(milo)
## ------------------------------------------------------------------
## 1) cell-level metadata 정리
##    - 기존 colData(milo)$sample 은 '환자 ID'라고 가정
##    - timepoint별 유일한 샘플 ID(sample_id)를 새로 만듦
## ------------------------------------------------------------------

colData(milo)$patient_id <- colData(milo)$sample
colData(milo)$sample_id  <- paste0(colData(milo)$patient_id, "_", colData(milo)$time)
meta$sample_id <- gsub('-','_',meta$sample_id)
## response 붙이기 (meta$sample 이 환자 ID라고 가정)
resp_map <- setNames(meta$`Binarized response`, meta$sample_id)
colData(milo)$response <- resp_map[as.character(colData(milo)$patient_id)]

meta_df <- as.data.frame(colData(milo))

## ------------------------------------------------------------------
## 2) countCells는 반드시 sample_id 기준으로
## ------------------------------------------------------------------

milo <- countCells(
  milo,
  meta.data = meta_df,
  sample = "sample_id"
)

dim(nhoodCounts(milo))
head(colnames(nhoodCounts(milo)))

## ------------------------------------------------------------------
## 3) sample-level design.df 만들기
##    - rownames는 sample_id
##    - patient_id는 paired blocking variable로 유지
## ------------------------------------------------------------------

design.df <- meta_df |>
  dplyr::select(sample_id, patient_id, time, response, dplyr::any_of("batch")) |>
  dplyr::distinct()

rownames(design.df) <- design.df$sample_id

## nhoodCounts 열 순서와 맞추기
design.df <- design.df[colnames(nhoodCounts(milo)), , drop = FALSE]

## factor 정리
design.df$sample_id  <- factor(design.df$sample_id)
design.df$patient_id <- factor(design.df$patient_id)
design.df$time       <- factor(design.df$time, levels = c("base", "1st"))
design.df$response   <- factor(design.df$response, levels = c("DCB", "NCB"))

stopifnot(all(rownames(design.df) == colnames(nhoodCounts(milo))))

## ------------------------------------------------------------------
## 4) 디자인 매트릭스
##    response main effect는 patient_id와 겹치므로
##    time:response 로 직접 두는 게 해석이 깔끔함
## ------------------------------------------------------------------
design_formula <- ~ patient_id + time 

## DCB 내 base vs 1st
design_dcb <- design.df |>
  dplyr::filter(response == "DCB") |>
  dplyr::group_by(patient_id) |>
  dplyr::filter(dplyr::n_distinct(time) == 2) |>
  dplyr::ungroup() |>
  as.data.frame()

keep_ids <- intersect(colnames(nhoodCounts(milo)), design_dcb$sample_id)

design_dcb <- design_dcb[match(keep_ids, design_dcb$sample_id), , drop = FALSE]
rownames(design_dcb) <- design_dcb$sample_id

## 핵심: unused factor levels 제거
design_dcb$patient_id <- droplevels(factor(design_dcb$patient_id))
design_dcb$time       <- droplevels(factor(design_dcb$time, levels = c("base", "1st")))
design_dcb$response   <- droplevels(factor(design_dcb$response))

milo_dcb <- milo
nhoodCounts(milo_dcb) <- nhoodCounts(milo)[, keep_ids, drop = FALSE]

stopifnot(identical(rownames(design_dcb), colnames(nhoodCounts(milo_dcb))))

mm_dcb <- model.matrix(~ patient_id + time, data = design_dcb)
colnames(mm_dcb)

da_dcb <- testNhoods(
  milo_dcb,
  design = ~ patient_id + time,
  design.df = design_dcb,
  model.contrasts = "time1st",
  subset.nhoods = NULL,
  fdr.weighting = "graph-overlap"
)

ggplot(da_dcb, aes(PValue))+geom_histogram()

milo_dcb <- buildNhoodGraph(milo_dcb)
da_dcb <- annotateNhoods(milo_dcb, da_dcb, coldata_col = "manual.cluster")
plotDAbeeswarm(da_dcb, group.by = "manual.cluster")
## NCB 내 base vs 1st
design_ncb <- design.df |>
  dplyr::filter(response == "NCB") |>
  dplyr::group_by(patient_id) |>
  dplyr::filter(dplyr::n_distinct(time) == 2) |>
  dplyr::ungroup() |>
  as.data.frame()

keep_ids <- intersect(colnames(nhoodCounts(milo)), design_ncb$sample_id)

design_ncb <- design_ncb[match(keep_ids, design_ncb$sample_id), , drop = FALSE]
rownames(design_ncb) <- design_ncb$sample_id

design_ncb$patient_id <- droplevels(factor(design_ncb$patient_id))
design_ncb$time       <- droplevels(factor(design_ncb$time, levels = c("base", "1st")))
design_ncb$response   <- droplevels(factor(design_ncb$response))

milo_ncb <- milo
nhoodCounts(milo_ncb) <- nhoodCounts(milo)[, keep_ids, drop = FALSE]

da_ncb <- testNhoods(
  milo_ncb,
  design = ~ patient_id + time,
  design.df = design_ncb,
  model.contrasts = "time1st",
  fdr.weighting = "graph-overlap"
)
ggplot(da_ncb, aes(PValue))+geom_histogram()

milo_ncb <- buildNhoodGraph(milo_ncb)
da_ncb <- annotateNhoods(milo_ncb, da_ncb, coldata_col = "manual.cluster")
plotDAbeeswarm(da_ncb, group.by = "manual.cluster")

## DCB+NCB 내 base vs 1st
stopifnot(identical(rownames(design.df), colnames(nhoodCounts(milo_dcb))))

mm <- model.matrix(~ time + (1|patient_id), data = design.df)
colnames(mm)

da_res <- testNhoods(
  milo,
  design = ~ time + (1|patient_id),
  design.df = design.df,
  model.contrasts = "time1st",
  subset.nhoods = NULL,
  fdr.weighting = "graph-overlap"
)
saveRDS(da_res, file='data/20260309_pilot/results/da_res.rds')
ggplot(da_res, aes(PValue))+geom_histogram()

milo <- buildNhoodGraph(milo)
da_res <- annotateNhoods(milo, da_res, coldata_col = "manual.cluster")
da_res <- annotateNhoods(milo, da_res, coldata_col = "manual.cluster_l2")
da_res$manual.cluster <- as.factor(da_res$manual.cluster)
plotDAbeeswarm(da_res, group.by = "manual.cluster")
head(da_res[!is.finite(da_res$logFC),])
bad_plot_idx <- !is.finite(da_res$logFC) | !is.finite(da_res$SpatialFDR)
da_res_clean <- da_res[!bad_plot_idx, , drop = FALSE]
plotDAbeeswarm(da_res_clean, group.by = "manual.cluster")

# ## 시간 변화가 DCB와 NCB 사이에서 다른가
# design.df <- as.data.frame(design.df)
# design.df$time <- factor(design.df$time, levels = c("base", "1st"))
# design.df$response <- factor(design.df$response, levels = c("DCB", "NCB"))
# design.df$patient_id <- factor(design.df$patient_id)
# 
# identical(rownames(design.df), colnames(nhoodCounts(milo)))
# 
# design_formula <- ~ time * response + (1 | patient_id)
# 
# ## 1) DCB 내 base vs 1st
# da_dcb_glmm <- testNhoods(
#   milo,
#   design = design_formula,
#   design.df = design.df,
#   model.contrasts = "time1st",
#   glmm.solver = "Fisher",
#   REML = TRUE,
#   fdr.weighting = "graph-overlap",
#   fail.on.error = FALSE
# )
# 
# ## 2) NCB 내 base vs 1st
# da_ncb_glmm <- testNhoods(
#   milo,
#   design = design_formula,
#   design.df = design.df,
#   model.contrasts = "time1st + time1st:responseNCB",
#   glmm.solver = "Fisher",
#   REML = TRUE,
#   fdr.weighting = "graph-overlap",
#   fail.on.error = FALSE
# )
# 
# ## 3) interaction: NCB의 변화량 - DCB의 변화량
# da_interaction_glmm <- testNhoods(
#   milo,
#   design = design_formula,
#   design.df = design.df,
#   model.contrasts = "time1st:responseNCB",
#   glmm.solver = "Fisher",
#   REML = TRUE,
#   fdr.weighting = "graph-overlap",
#   fail.on.error = FALSE
# )
# 
# save(da_dcb, da_ncb, da_dcb_glmm, file = "data/20260309_pilot/results/20260408_milo_results.RData")