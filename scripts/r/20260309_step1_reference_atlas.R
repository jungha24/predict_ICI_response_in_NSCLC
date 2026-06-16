### 20260309
### make annotation reference atlas
### steps required fine-tuning "#<<<"

### install
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))
remotes::install_github("bnprks/BPCells/r")
remotes::install_version(
  package = "SeuratObject",
  version = "5.0.0",
  repos = "https://cloud.r-project.org"
)
devtools::install_github("satijalab/seurat-data", "seurat5")
#devtools::install_github("satijalab/azimuth", "seurat5")
remotes::install_github('satijalab/azimuth', ref = 'master')

# aligning version of Seurat and SeuratObject
remove.packages("SeuratObject")
install.packages("SeuratObject", repos = "https://cloud.r-project.org")
remove.packages("Seurat")
install.packages("Seurat", repos = "https://cloud.r-project.org")

### package
library(reticulate)
use_python("/opt/venvs/seurat-py/bin/python", required = TRUE)
library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
#library(SeuratDisk)
library(BPCells)
library(SeuratObject)
library(Matrix)
library(sctransform)
library(Azimuth)
library(SeuratData)
library(SeuratObject)
library(patchwork)
library(future)
library(ggplot2)
library(stringr)
library(ggthemes)
library(dplyr)
library(ggbeeswarm)
library(reshape)
library(tidyr)
library(arrow)
library(Seurat)
library(tibble)

### environment
plan("sequential")
options(future.globals.maxSize = Inf)

### [1] load data
## 1. nsclc data
sce <- readH5AD('data/20260309_pilot/nsclc_n73/20260309_adata_QC.h5ad')
assayNames(sce)
obj <- CreateSeuratObject(
  counts = assay(sce,"counts"),
  meta.data = as.data.frame(colData(sce))
)

obj$cohort <- 'Lee_p1'
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mt",'pct_counts_ribo','pct_counts_hb'), ncol = 5,alpha=0)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mt",'pct_counts_ribo','pct_counts_hb'), ncol = 5,alpha=0)
# RBC check
mat <- LayerData(obj, assay = "RNA", layer = "counts")
hb_sum <- mat["HBA1", ] + mat["HBB", ]
sum(hb_sum>10)
obj <- subset(obj, cells = colnames(obj)[hb_sum <= 10])
saveRDS(obj,file='data/20260309_pilot/nsclc_n73/20260309_adata_QC.rds')

# ## 2. covid control data
# cov_hc <- readRDS('data/250729_COVID_final.rds') # covid normal
# cov_hc <- subset(aida, subset = severity =="CTRL") # covid normal
# raw_counts_normal <- LayerData(cov_hc[["RNA"]], layer = "counts.1")
# cov_hc_raw <- CreateSeuratObject(
#   counts = raw_counts_normal,
#   project = "HC_covidCohort"
# )

## 2. AIDA control data
mat <- readMM("data/20260309_pilot/aida/aida_korean_matrix.mtx") # generated from data/20260309_pilot/test.ipynb
obs <- read.csv("data/20260309_pilot/aida/aida_korean_obs.csv", row.names = 1) # generated from data/20260309_pilot/test.ipynb
var <- read.csv("data/20260309_pilot/aida/aida_korean_var.csv", row.names = 1) # generated from data/20260309_pilot/test.ipynb
# AnnData는 X가 cells x genes 이므로, Seurat용으로 전치
mat <- t(mat)
rownames(mat) <- rownames(var)
colnames(mat) <- rownames(obs)
aida <- CreateSeuratObject(
  counts = mat,
  meta.data = obs
)
aida@misc$feature_metadata <- var
aida
saveRDS(aida, file='data/20260309_pilot/aida/aida_korean.rds')

# feature name이 ensembl로 되어 있었음
aida@misc$ensembl_ids <- rownames(aida)
gene_symbols <- aida@misc$feature_metadata$feature_name
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- rownames(aida)[is.na(gene_symbols) | gene_symbols == ""]
gene_symbols <- make.unique(gene_symbols)
rownames(aida) <- gene_symbols

aida[["pct_counts_mt"]] <- PercentageFeatureSet(aida, pattern = "^MT-")
aida[["pct_counts_ribo"]] <- PercentageFeatureSet(aida, pattern = "^RP[SL]")
aida[["pct_counts_hb"]]   <- PercentageFeatureSet(aida, pattern = "^HB[AB]")
VlnPlot(aida, features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mt",'pct_counts_ribo','pct_counts_hb'), ncol = 5,alpha=0) #3x8 #20260310_aida_korean_QC
# RBC check
mat <- LayerData(aida, assay = "RNA", layer = "counts")
hb_sum <- mat["HBA1", ] + mat["HBB", ]
sum(hb_sum>10)
aida <- subset(aida, subset = nCount_RNA > 1800 & nCount_RNA < 8000 &  nFeature_RNA < 4000 & pct_counts_mt < 2.5 )
VlnPlot(aida, features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mt",'pct_counts_ribo','pct_counts_hb'), ncol = 5,alpha=0) #3x8 #20260310_aida_korean_QC_mask
aida$cohort <- 'AIDA'
saveRDS(aida,file='data/20260309_pilot/aida/20260310_aida_QC.rds')

### [2] LogNorm
obj <- RenameCells(obj, add.cell.id = "DIS")
aida  <- RenameCells(aida,  add.cell.id = "NOR")

obj_list <- list(
  disease = obj,
  normal  = aida
)
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) #<<<
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000) #<<<
  return(x)
})

seu_integrated$cohort_l2 <- paste0(seu_integrated$time, '_',seu_integrated$cohort)
unique(seu_integrated@meta.data[seu_integrated@meta.data$time =='1st',]$sample)
length(unique(seu_integrated@meta.data[seu_integrated@meta.data$cohort =='Lee_p1',]$sample))
length(unique(seu_integrated@meta.data[seu_integrated@meta.data$cohort =='Lee_p1' & seu_integrated$time=='base',]$sample))
length(unique(seu_integrated@meta.data[seu_integrated@meta.data$cohort =='Lee_p1' & seu_integrated$time=='1st',]$sample))
length(unique(seu_integrated@meta.data[seu_integrated@meta.data$cohort =='Lee_p1' & seu_integrated$time=='2nd',]$sample))
DimPlot(seu_integrated, reduction = "umap", group.by='seurat_clusters') #4x5 #20260311_integrated_res0_15 #20260311_integrated_pca20_res0_2 #20260311_integrated_pca20_res0_3
p <- DimPlot(seu_integrated, reduction = "umap", group.by='integrated_snn_res.0.2', split.by = 'integrated_snn_res.0.2')
p + facet_wrap(~integrated_snn_res.0.2)#12x14 #20260311_integrated_res0_15_splitBySeuratcluster #20260311_integrated_pca20_res0_2_splitBySeuratcluster
DimPlot(seu_integrated, reduction = "umap", group.by='seurat_clusters', split.by='cohort_l2')#20260311_integrate_pca20_res03_splitBycohort
p <- DimPlot(seu_integrated, reduction = "umap", group.by='seurat_clusters', split.by='seurat_clusters')
p+facet_wrap(~seurat_clusters)#12x14 #20260311_integrate_pca20_res03_splitBySeuratcluster

### [5] cell type annotation
### [5-1] azimuth annotation
DefaultAssay(seu_integrated) <- "RNA"
seu_integrated[["percent.mt"]] <- PercentageFeatureSet(seu_integrated, pattern = "^MT-")
#seu_integrated[["RNA"]] <- JoinLayers(seu_integrated) # wrong method
seu_integrated <- JoinLayers(seu_integrated, assay = "RNA")
Layers(seu_integrated[["RNA"]])
seu_integrated <- RunAzimuth(seu_integrated, reference = "pbmcref")

DefaultAssay(seu_integrated) <- 'integrated'
DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l1', split.by = 'cohort_l2')#4x11 #20260311_integrated_azimuth_l1_splitByCohort #20260311_integrated_pca20_azimuth_l1_splitByCohort
DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l2', split.by = 'cohort_l2')#4x11 #20260311_integrated_azimuth_l2_splitByCohort #20260311_integrated_pca20_azimuth_l2_splitByCohort
p<- DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l2',split.by='predicted.celltype.l2')
p + facet_wrap(~predicted.celltype.l2)#12x14 #20260311_integrated_azimuth_l2 #20260311_integrated_pca20_azimuth_l2
saveRDS(seu_integrated, file='data/20260309_pilot/results/20260310_integrated.rds')

### [5-2] CITE-seq reference containing 162,000 cells and 228 antibodies (https://zenodo.org/records/7779017/files/pbmc_multimodal_2023.rds)
cite_ref <- readRDS('data/20260309_pilot/pbmc_multimodal_2023.rds') # azimuth와 같은 데이터인데 좀 더 latest느낌남.
DefaultAssay(seu_integrated) <- "RNA"
anchor <- FindTransferAnchors(
  reference = cite_ref,
  query = seu_integrated,
  reference.reduction='spca',
  normalization.method = 'SCT',
  dims=1:50
)
seu_integrated <- MapQuery(
  anchorset = anchor,
  query = seu_integrated,
  reference = cite_ref,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"
  ),
  reduction.model = "wnn.umap"
)

# azimuth의 predicted.celltype.l1.score/ predicted.celltype.l1/predicted.celltype.l2.score/predicted.celltype.l2가 덮어써짐...
DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l1', split.by = 'cohort_l2')#4x11 #20260311_integrated_pca20_cite_l1_splitByCohort
DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l2', split.by = 'cohort_l2')#4x11 #20260311_integrated_pca20_cite_l2_splitByCohort
p<- DimPlot(seu_integrated, reduction = "umap", group.by='predicted.celltype.l2',split.by='predicted.celltype.l2')
p + facet_wrap(~predicted.celltype.l2)#12x14 #20260311_integrated_pca20_cite_l2

### [5-3] manual curation
# integrated assay: PCA PC 1-20; res 0.2

cells_use <- rownames(seu_integrated@meta.data)[
  seu_integrated$file == "9-IO_1st"
]
DimPlot(seu_integrated,
        reduction = "umap",
        cells.highlight = cells_use,
        cols.highlight = "red",
        cols="grey80", order=T)


seu_integrated$time <- str_split_fixed(seu_integrated$file,'_',2)[,2]
seu_integrated$cohort_l2 <- paste0(seu_integrated$cohort, '_', seu_integrated$time)

ggplot(seu_integrated@meta.data, aes(x=seurat_clusters, fill=cohort_l2))+geom_bar()+theme_bw()#3x5#20260312_integrated_pca20_res02_bar_byCohort
ggplot(seu_integrated@meta.data, aes(x=seurat_clusters, fill=cohort_l2))+geom_bar(position='fill')+theme_bw()#3x5#20260312_integrated_pca20_res02_bar_byCohort_fill

DimPlot(seu_integrated, reduction = "umap", group.by='seurat_clusters',label=T)#3x4 #20260312_integrated_pca20_res02_label
p <- DimPlot(seu_integrated, reduction = "umap", group.by='seurat_clusters', split.by='seurat_clusters')
p+facet_wrap(~seurat_clusters)#12x14 #20260312_integrated_pca20_res02_splitBySeuratclusters
tab_cohort <- table(seu_integrated$seurat_clusters, seu_integrated$cohort_l2)
prop_cohort_by_cluster <- prop.table(tab_cohort, margin = 1)
round(prop_cohort_by_cluster, 3)

markers <- FindAllMarkers(seu_integrated, only.pos = TRUE)
markers %>%
  group_by(seurat_clusters) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(seu_integrated, features = top10$gene) + NoLegend() #13x18 #20260312_integrated_pca20_res02_top10marker_heatmap

markers.pick <- c('TCF7', 'LEF1', 'CCR7', 'IL7R','KLRG1',
  'NKG7', 'CCL5', 'KLRD1', 'KLRF1','FGFBP2', 'PRF1', 'CTSW', 'GNLY','GZMB', 'GZMH',
  'CD8A', 'CD8B', 'CD3D', 'CD3G',
  'CD3D', 'CD3G', 'SIT1',
  'FCGR3A','S100A8', 'CTSS','LST1', 'CD14', 'AIF1', 'GRN', 'SPI1', 'CEBPD','CDKN1C',
  'CSF1R', 'MS4A7', 'LILRB2', 'MAFB','CDKN1C',
  'FCER1A', 'CD1C', 'CLEC10A', 'FLT3',
  'MS4A1', 'CD79A', 'BANK1', 'FCRL1', 'SPIB',
  'PPBP', 'PF4', 'TUBB1', 'MPIG6B', 'GP9')

seu_integrated$manual.cluster <-
  ifelse(seu_integrated$seurat_clusters %in% c(0, 14),'CD14 Mono',
         ifelse(seu_integrated$seurat_clusters == 1, 'NK',
                ifelse(seu_integrated$seurat_clusters == 2,'effector CD8+T',
                       ifelse(seu_integrated$seurat_clusters %in% c(3,4), 'naive/TCM CD4+T',
                              ifelse(seu_integrated$seurat_clusters %in% c(5,12), 'B',
                                     ifelse(seu_integrated$seurat_clusters ==6, 'CD16 Mono',
                                            ifelse(seu_integrated$seurat_clusters ==7, 'navie CD8+T',
                                                   ifelse(seu_integrated$seurat_clusters %in% c(8,10), 'platelet',
                                                          ifelse(seu_integrated$seurat_clusters ==9,'T/NK',
                                                                 ifelse(seu_integrated$seurat_clusters == 11,'cDC2','ambiguous'))))))))))

DotPlot(seu_integrated, features = unique(markers.pick), group.by= 'manual.cluster', cols=c('white','red')) + RotatedAxis()#4x12 ##20260312_integrated_pca20_manualcluster_dot_marker
Idents(seu_integrated) <- 'manual.cluster'
DoHeatmap(seu_integrated, features = unique(markers.pick)) + NoLegend() #13x18 #20260312_integrated_pca20_manualcluster_heatmap
DimPlot(seu_integrated,reduction = "umap",  group.by='manual.cluster',label=T) #20260312_integrated_pca20_manualcluster_umap
DimPlot(seu_integrated,reduction = "umap",  group.by='manual.cluster',split.by = 'cohort_l2',label=T) #20260312_integrated_pca20_manualcluster_splitByCohortl2
ggplot(seu_integrated@meta.data, aes(x=manual.cluster,fill=cohort_l2))+geom_bar()+theme_bw()+theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5)) #4x5 #20260312_integrated_pca20_manualcluster_byCohortl2
ggplot(seu_integrated@meta.data, aes(x=manual.cluster,fill=cohort_l2))+geom_bar(position='fill')+theme_bw()+theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5)) #4x5 #20260312_integrated_pca20_manualcluster_byCohortl2_v2

 
df <- as.data.frame(table(seu_integrated$manual.cluster, seu_integrated$sample_time))
df$donor <- str_split_fixed(df$Var2,'-',2)[,1]
df$time <- str_split_fixed(df$Var2,'-',2)[,2]
df$time <- ifelse(df$time =='2nd','1st',df$time)
ggplot(df, aes(x=time, y=Freq))+geom_boxplot(outlier.size = 0)+geom_beeswarm(size=0.6,alpha=0.3,cex=2)+facet_wrap(~Var1, scales='free_y')+theme_few() #6x8 #20260312_integrated_pca20_manualcluster_bySample

### [6] robustness 확인
### 1. 환자별 cell type count/ proportion 계산
### nsclc_n73 all, base only, post only 나눠서 계산
cell_df <- seu_integrated@meta.data[seu_integrated@meta.data$time =='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  dplyr::select(sample, manual.cluster)
count_df <- seu_integrated@meta.data[seu_integrated@meta.data$time =='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  count(sample, manual.cluster, name = "n_cells")
total_df <- seu_integrated@meta.data[seu_integrated@meta.data$time =='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  count(sample, name = "total_cells")
prop_df <- count_df %>%
  left_join(total_df, by = "sample") %>%
  mutate(prop = n_cells / total_cells)

cell_df <- seu_integrated@meta.data[seu_integrated@meta.data$time !='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  dplyr::select(sample, manual.cluster)
count_df <- seu_integrated@meta.data[seu_integrated@meta.data$time !='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  count(sample, manual.cluster, name = "n_cells")
total_df <- seu_integrated@meta.data[seu_integrated@meta.data$time !='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  count(sample, name = "total_cells")
prop_df <- count_df %>%
  left_join(total_df, by = "sample") %>%
  mutate(prop = n_cells / total_cells)

### 2. prevalence/median proportion/CV/ 최소 cell 수
gini <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (sum(x) == 0) return(NA_real_)
  (2 * sum(x * seq_len(n)) / (n * sum(x))) - (n + 1) / n
}

robust_df <- prop_df %>%
  group_by(manual.cluster) %>%
  summarise(
    n_patients_observed = n_distinct(sample),
    prevalence = n_distinct(sample) / n_distinct(cell_df$sample),
    median_prop = median(prop, na.rm = TRUE),
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    cv_prop = sd_prop / mean_prop,
    gini_prop = gini(prop),
    median_n_cells = median(n_cells, na.rm = TRUE),
    pct_patients_ge10 = mean(n_cells >= 10),
    pct_patients_ge20 = mean(n_cells >= 20),
    pct_patients_ge50 = mean(n_cells >= 50),
    pct_patients_ge100 = mean(n_cells >= 100),
  ) %>%
  arrange(desc(prevalence), desc(median_n_cells))

robust_df
## zero 포함해서 더 정확하게
all_patients <- unique(cell_df$sample)
all_celltypes <- unique(cell_df$manual.cluster)

full_df <- tidyr::expand_grid(sample = all_patients, manual.cluster = all_celltypes) %>%
  left_join(count_df, by = c("sample", "manual.cluster")) %>%
  left_join(total_df, by = "sample") %>%
  mutate(
    n_cells = ifelse(is.na(n_cells), 0, n_cells),
    prop = n_cells / total_cells
  )

robust_df2 <- full_df %>%
  group_by(manual.cluster) %>%
  summarise(
    prevalence = mean(n_cells > 0),
    zero_fraction = mean(n_cells == 0),
    median_prop = median(prop, na.rm = TRUE),
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    cv_prop = ifelse(mean_prop == 0, NA, sd_prop / mean_prop),
    gini_prop = gini(prop),
    median_n_cells = median(n_cells, na.rm = TRUE),
    mean_n_cells = median(n_cells, na.rm = TRUE),
    minimum_n_cells = min(n_cells, na.rm = TRUE),
    max_n_cells = max(n_cells, na.rm = TRUE),
    sd_n_cells = sd(n_cells, na.rm=TRUE),
    cv_n_cells = ifelse(mean_n_cells == 0, NA, sd_n_cells / mean_n_cells),
    gini_n_cells = gini(n_cells),
    pct_patients_ge10 = mean(n_cells >= 10),
    pct_patients_ge20 = mean(n_cells >= 20),
    pct_patients_ge50 = mean(n_cells >= 50),
    pct_patients_ge100 = mean(n_cells >= 100)
  ) %>%
  arrange(desc(prevalence), desc(median_n_cells))

robust_df2$time <- 'base'
robust_df3$time <- 'post'

robust_df4 <- rbind(robust_df2,robust_df3)
write.table(robust_df4, file='data/20260309_pilot/results/20260312_integrated_pca20_manual.cluster_celltype_spec_table.txt',quote=F, sep='\t')
robust_df4 <- as.data.frame(robust_df4)

robust_df4.melt <- robust_df4 |>
  pivot_longer(
    cols = -c(manual.cluster, time),
    names_to = "variable",
    values_to = "value"
  )
robust_df4.melt$manual.cluster <- factor(robust_df4.melt$manual.cluster,
                                         levels=c('CD14 Mono','naive/TCM CD4+T','NK','effector CD8+T','B','CD16 Mono','T/NK','navie CD8+T','cDC2','platelet','ambiguous'))
robust_df4.melt$variable <- factor(robust_df4.melt$variable, 
                                   levels=c('prevalence','zero_fraction','median_prop','mean_prop','sd_prop','cv_prop','gini_prop',
                                            'median_n_cells','mean_n_cells','sd_n_cells','cv_n_cells','gini_n_cells','minimum_n_cells','max_n_cells',
                                            'pct_patients_ge10','pct_patients_ge20','pct_patients_ge50','pct_patients_ge100'))
 
ggplot(robust_df4.melt, aes(x=manual.cluster, y=value,col=time,group=time))+geom_line()+facet_wrap(~variable, scales='free_y',ncol=7)+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5)) #6x8 #20260312_integrated_pca20_manual.cluster_celltype_spec
save(seu_integrated, file='data/20260309_pilot/results/20260310_integrated.RData')

ggplot(full_df,aes(x='celltype',y=n_cells))+geom_boxplot(outlier.size = 0)+geom_beeswarm(cex=3, size=0.6, alpha=0.4)+facet_wrap(.~factor(manual.cluster,levels=c('CD14 Mono','naive/TCM CD4+T','NK','effector CD8+T','B','CD16 Mono','T/NK','navie CD8+T','cDC2','platelet','ambiguous')), scales='free_y')+ylim(c(0,30))#6x6 #20260313_integrated_pca20_manualcluster_n_cells

patients_ge10_by_celltype <- count_df %>%
  filter(n_cells >= 10) %>%
  group_by(manual.cluster) %>%
  summarise(patient_set = list(sort(unique(sample))), .groups='drop')

patient_sets_list <- setNames(
  patients_ge10_by_celltype$patient_set,
  patients_ge10_by_celltype$manual.cluster
)
selected_celltypes <- c('B', 'CD14 Mono', 'CD16 Mono', 'effector CD8+T', 'naive/TCM CD4+T', 'NK', 'T/NK')
selected_sets <- patient_sets_list[selected_celltypes]
common_patients_selected <- Reduce(intersect, selected_sets)
length(common_patients_selected)

df_common <- count_df %>%
  filter(sample %in% common_patients_selected,
         manual.cluster %in% selected_celltypes)
df_common %>%
  arrange(manual.cluster, n_cells)

min_ncells_by_celltype <- df_common %>%
  group_by(manual.cluster) %>%
  summarise(
    min_n_cells = min(n_cells),
    median_n_cells = median(n_cells),
    max_n_cells = max(n_cells),
    .groups = "drop"
  ) %>%
  arrange(min_n_cells)

min_ncells_by_celltype

patients_le20_by_celltype <- df_common %>%
  filter(n_cells <= 20) %>%
  group_by(manual.cluster) %>%
  summarise(
    n_patients_le20 = n(),
    patient_ids = paste(sort(sample), collapse = ", "),
    .groups = "drop"
  )

patients_le20_by_celltype

#common_patients_selected[!(common_patients_selected %in% c('IO_SC_021', 'IO_SC_037', 'IO_SC_064', 'IO_SC_068','IO_SC_033', 'IO_SC_080', 'IO_SC_126'))]

# 결국 B, CD16 Mono는 PCA안쓰기로 metadata에 없는 sample이 있어걔들빼면 몇안남을거같아 
length(common_patients_selected[!(common_patients_selected %in% c('IO_SC_005','IO_SC_032','IO_SC_041','IO_SC_053','IO_SC_054','IO_SC_061','IO_SC_081','IO_SC_089','IO_SC_091','IO_SC_097','IO_SC_104','IO_SC_110','Normal'))])

common_patients_selected <- common_patients_selected[!(common_patients_selected %in% c('IO_SC_005','IO_SC_032','IO_SC_041','IO_SC_053','IO_SC_054','IO_SC_061','IO_SC_081','IO_SC_089','IO_SC_091','IO_SC_097','IO_SC_104','IO_SC_110','Normal'))]
write.table(common_patients_selected, file='data/20260309_pilot/results/20260313_patients_selected.txt',row.names = F, col.names=F, quote=F)

### [6] integrated annotation을 raw object에 옮기기 
obj <- readRDS('data/20260309_pilot/nsclc_n73/20260309_adata_QC.rds')
aida <- readRDS('data/20260309_pilot/aida/20260310_aida_QC.rds')

obj <- RenameCells(obj, add.cell.id = "DIS")
aida  <- RenameCells(aida,  add.cell.id = "NOR")
aida$sample <- aida$donor_id
aida$Age_IO <- as.numeric(str_split_fixed(aida$development_stage,'-',3)[,1])
aida$Sex.M..1..F..2 <- ifelse(aida$sex =='male',1,2)
raw_query <- merge(
  x = obj,
  y = aida,
  project = "PBMC_raw_merged"
)

common_cells <- intersect(colnames(raw_query), colnames(seu_integrated))
length(common_cells)

raw_query$celltype_main <- NA_character_
raw_query$celltype_main[match(common_cells, colnames(raw_query))] <-
  seu_integrated[['manual.cluster']][match(common_cells, colnames(seu_integrated)), 1]

table(is.na(raw_query$celltype_main))
table(raw_query$celltype_main, useNA = "ifany")

saveRDS(raw_query, file='data/20260309_pilot/results/20260313_raw_query_annotation_attached.rds')

### [7] Discovery (step2)에 필요한 데이터 저장하기
raw_query <- JoinLayers(raw_query, assay = "RNA")
counts <- GetAssayData(raw_query, assay = "RNA", layer = "counts")

# MatrixMarket export
writeMM(counts, "data/20260309_pilot/results/version1/data_export/20260313_counts.mtx")
write.table(
  rownames(counts),
  file = "data/20260309_pilot/results/version1/data_export/20260313_genes.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

write.table(
  colnames(counts),
  file = "data/20260309_pilot/results/version1/data_export/20260313_barcodes.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

cell_meta <- raw_query@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(
    patient_id = as.character(sample),
    source = as.character(cohort),
    collection_batch = as.character(file),
    timepoint = as.character(time),
    celltype_main = as.character(celltype_main),
    nCount_RNA = as.numeric(nCount_RNA),
    Sample.ID = as.character(Sample.ID),
    Study.ID = as.character(Study.ID),
    Age_IO = as.numeric(Age_IO),
    sex = as.numeric(Sex.M..1..F..2),
    Histology = as.character(Histology),
    smoking = as.character(Smoking.Never..0.Ex..1.Current..2),
    ecog = as.numeric(ECOG._PS),
    drug = as.character(Drug),
    PD_event = as.numeric(PD_Event),
    PFS = as.numeric(PFS..Days.),
    Death_Event = as.numeric(Death_Event),
    OS = as.numeric(OS..Days.),
    EGFR_Mutation_Status = as.character(EGFR_Mutation_Status),
    IO_Line = as.character(IO_Line),
    Previous_palliative_chemo = as.character(Previous_palliative_chemo),
    Previous_palliative_target = as.character(Previous_palliative_target),
    PD.L1_TPS = as.character(PD.L1_TPS)
    
  ) %>%
  select(cell_id, patient_id, source, collection_batch, timepoint, celltype_main, nCount_RNA, Sample.ID, Study.ID, Age_IO, sex,Histology, smoking,ecog,drug, PD_event, PFS, Death_Event, OS, EGFR_Mutation_Status, IO_Line, Previous_palliative_chemo, Previous_palliative_target, PD.L1_TPS)

arrow::write_parquet(cell_meta, "data/20260309_pilot/results/version1/data_export/20260313_cell_meta.parquet")
### save sessionInfo
ip <- installed.packages()
pkg_df <- data.frame(
  Package = ip[, "Package"],
  Version = ip[, "Version"],
  LibPath = ip[, "LibPath"],
  row.names = NULL
)

write.csv(pkg_df, "data/20260309_pilot/installed_packages_versions_step1.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "data/20260309_pilot/sessionInfo_step1.txt")
