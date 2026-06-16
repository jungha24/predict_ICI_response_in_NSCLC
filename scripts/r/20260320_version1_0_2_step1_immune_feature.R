# 2026.03.20
# make biological meaningful SBP

load('data/20260309_pilot/results/20260310_integrated.RData')

# myleid: CD14 Mono, CD16 Mono, cDC2
# lymphoid: B, effector CD8+T, navie/TCM CD4+T, navie CD8+T, NK, T/NK

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tibble)

# 1. parameters
sbp_pseudocount <- 0.5

myeloid_types <- c(
  "CD14 Mono",'CD16 Mono','cDC2'
)
lymphoid_types <- c(
  'B','effector CD8+T','naive/TCM CD4+T','navie CD8+T','NK','T/NK'
)

mono_balance_types <- c("CD14 Mono", "CD16 Mono")

naive_genes <- c("TCF7", "LEF1", "CCR7", "IL7R", "MAL", "LTB")
cytotoxic_genes <- c("NKG7", "CCL5", "KLRD1", "KLRF1", "FGFBP2", "PRF1", "CTSW", "GNLY", "GZMB", "GZMH", "CCL4")
cd14_genes <- c("CD14", "VCAN", "S100A8", "S100A9", "LYZ")
cd16_genes <- c("FCGR3A", "LST1", "CDKN1C", "MS4A7", "IFITM3", "CX3CR1", "GPBAR1", "LRRC25", "CALHM6", "CSF1R", "MS4A7", "LILRB2", "MAFB")

selected_patients <- fread('data/20260309_pilot/results/20260313_patients_selected.txt',header=F)

# 2. patient-level L/M balance
lm_count_long <- seu_integrated@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  filter(sample %in% selected_patients$V1, time == "base") %>%
  select(sample, manual.cluster) %>%
  mutate(
    lm_group = case_when(
      manual.cluster %in% lymphoid_types ~ "lymphoid",
      manual.cluster %in% myeloid_types  ~ "myeloid",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(lm_group)) %>%
  count(sample, lm_group, name = "n_cells")

lm_patient_df <- lm_count_long %>%
  tidyr::pivot_wider(
    names_from = lm_group,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  mutate(
    lymphoid = ifelse(is.na(lymphoid), 0L, lymphoid),
    myeloid  = ifelse(is.na(myeloid),  0L, myeloid),
    lm_total = lymphoid + myeloid,
    lymphoid_prop = ifelse(lm_total > 0, lymphoid / lm_total, NA_real_),
    myeloid_prop  = ifelse(lm_total > 0, myeloid  / lm_total, NA_real_),
    LM_log_ratio = log((lymphoid + sbp_pseudocount) / (myeloid + sbp_pseudocount)),
    LM_SBP_2part_ilr = sqrt(1/2) * LM_log_ratio,
    LM_low_total_flag = lm_total < 50
  ) %>%
  rename(patient_id = sample) %>%
  arrange(patient_id)


# patient-level table도 따로 저장
write.table(
  lm_patient_df,
  file = "data/20260309_pilot/results/version1/data_export/20260320_patient_LM_SBP.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

# 3. patient-level CD14/DC16 Mono balance
mono_count_long <- seu_integrated@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  filter(sample %in% selected_patients$V1, time == "base") %>%
  filter(manual.cluster %in% mono_balance_types) %>%
  count(sample, manual.cluster, name = "n_cells") %>%
  mutate(
    mono_group = case_when(
      manual.cluster == "CD14 Mono" ~ "cd14_mono",
      manual.cluster == "CD16 Mono" ~ "cd16_mono",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(mono_group)) %>%
  select(sample, mono_group, n_cells)

mono_patient_df <- mono_count_long %>%
  tidyr::pivot_wider(
    names_from = mono_group,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  mutate(
    cd14_mono = ifelse(is.na(cd14_mono), 0L, cd14_mono),
    cd16_mono = ifelse(is.na(cd16_mono), 0L, cd16_mono),
    mono_total = cd14_mono + cd16_mono,
    cd14_mono_prop = ifelse(mono_total > 0, cd14_mono / mono_total, NA_real_),
    cd16_mono_prop = ifelse(mono_total > 0, cd16_mono / mono_total, NA_real_),
    CD14_CD16_log_ratio = log((cd14_mono + sbp_pseudocount) / (cd16_mono + sbp_pseudocount)),
    CD14_CD16_SBP_2part_ilr = sqrt(1/2) * CD14_CD16_log_ratio,
    CD14_CD16_low_total_flag = mono_total < 20
  ) %>%
  rename(patient_id = sample) %>%
  arrange(patient_id)

write.table(
  mono_patient_df,
  file = "data/20260309_pilot/results/version1/data_export/20260320_patient_CD14_CD16_SBP.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

# 4.raw_query RNA assay에서 addmodulescore

score_obj <- subset(obj, subset = sample %in% selected_patients$V1 & time =='base')
DefaultAssay(score_obj) <- "RNA"
score_obj <- NormalizeData(score_obj, verbose = FALSE)
score_obj <- RenameCells(score_obj, add.cell.id = "DIS")
common_cells <- intersect(colnames(score_obj), colnames(seu_integrated))
length(common_cells)

score_obj$celltype_main <- NA_character_
score_obj$celltype_main[match(common_cells, colnames(score_obj))] <-
  seu_integrated[['manual.cluster']][match(common_cells, colnames(seu_integrated)), 1]

# 5. helper function
safe_gene_set <- function(obj, genes, label = "gene_set") {
  genes_found <- intersect(unique(genes), rownames(obj))
  if (length(genes_found) == 0) {
    stop(paste0("No genes from ", label, " found in object."))
  }
  if (length(genes_found) < length(unique(genes))) {
    missing_genes <- setdiff(unique(genes), genes_found)
    message(label, " missing genes dropped: ", paste(missing_genes, collapse = ", "))
  }
  genes_found
}

add_subset_module_scores <- function(obj, subset_clusters, feature_list, seed = 1) {
  meta0 <- obj@meta.data %>%
    tibble::rownames_to_column("cell_id")
  
  subset_cells <- meta0 %>%
    filter(celltype_main %in% subset_clusters) %>%
    pull(cell_id)
  
  sub_obj <- subset(obj, cells = subset_cells)
  
  for (nm in names(feature_list)) {
    sub_obj <- AddModuleScore(
      object = sub_obj,
      features = list(feature_list[[nm]]),
      name = paste0(nm, "_"),
      assay = "RNA",
      search = FALSE,
      seed = seed
    )
    score_col <- paste0(nm, "_1")
    sub_obj@meta.data[[nm]] <- sub_obj@meta.data[[score_col]]
    sub_obj@meta.data[[score_col]] <- NULL
  }
  
  score_df <- sub_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    select(cell_id, all_of(names(feature_list)))
  
  meta0 %>%
    select(cell_id) %>%
    left_join(score_df, by = "cell_id")
}

# 6. gene set 정리
naive_genes_use <- safe_gene_set(raw_query, naive_genes, "naive_genes")
cytotoxic_genes_use <- safe_gene_set(raw_query, cytotoxic_genes, "cytotoxic_genes")
cd14_genes_use <- safe_gene_set(raw_query, cd14_genes, "cd14_genes")
cd16_genes_use <- safe_gene_set(raw_query, cd16_genes, "cd16_genes")

# 7. lymphoid subset에서 naive / cytotoxic score 계산
t_lymhpoid_types = c("effector CD8+T","naive/TCM CD4+T","navie CD8+T","NK","T/NK")
lymphoid_score_df <- add_subset_module_scores(
  obj = score_obj,
  subset_clusters = t_lymhpoid_types,
  feature_list = list(
    naive_score = naive_genes_use,
    cytotoxic_score = cytotoxic_genes_use
  ),
  seed = 1
) %>%
  mutate(
    naive_minus_cytotoxic_score = naive_score - cytotoxic_score
  )

# 8. CD14/CD16 monocyte subset에서 score 계산
mono_score_df <- add_subset_module_scores(
  obj = score_obj,
  subset_clusters = mono_balance_types,
  feature_list = list(
    cd14_score = cd14_genes_use,
    cd16_score = cd16_genes_use
  ),
  seed = 1
) %>%
  mutate(
    cd14_minus_cd16_score = cd14_score - cd16_score
  )

# 9. cell-level score table 저장
cell_score_df <- score_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id) %>%
  left_join(
    lymphoid_score_df %>%
      select(cell_id, naive_score, cytotoxic_score, naive_minus_cytotoxic_score),
    by = "cell_id"
  ) %>%
  left_join(
    mono_score_df %>%
      select(cell_id, cd14_score, cd16_score, cd14_minus_cd16_score),
    by = "cell_id"
  )

write.table(
  cell_score_df,
  file = "data/20260309_pilot/results/version1/data_export/20260320_cell_level_AddModuleScores.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

# 10. patient-level score summary (Lee_p1 + base only)
cell_to_patient_df <- score_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  transmute(
    cell_id,
    patient_id = as.character(sample),
    cohort = as.character(cohort),
    time = as.character(time)
  )

lymphoid_score_patient_df <- lymphoid_score_df %>%
  left_join(cell_to_patient_df, by = "cell_id") %>%
  filter(!is.na(naive_score) | !is.na(cytotoxic_score)) %>%
  group_by(patient_id) %>%
  summarise(
    patient_n_lymphoid_scored_cells = n(),
    patient_naive_score_mean = mean(naive_score, na.rm = TRUE),
    patient_cytotoxic_score_mean = mean(cytotoxic_score, na.rm = TRUE),
    patient_naive_minus_cytotoxic_score_mean = mean(naive_minus_cytotoxic_score, na.rm = TRUE),
    .groups = "drop"
  )


mono_score_patient_df <- mono_score_df %>%
  left_join(cell_to_patient_df, by = "cell_id") %>%
  filter(!is.na(cd14_score) | !is.na(cd16_score)) %>%
  group_by(patient_id) %>%
  summarise(
    patient_n_mono_scored_cells = n(),
    patient_cd14_score_mean = mean(cd14_score, na.rm = TRUE),
    patient_cd16_score_mean = mean(cd16_score, na.rm = TRUE),
    patient_cd14_minus_cd16_score_mean = mean(cd14_minus_cd16_score, na.rm = TRUE),
    .groups = "drop"
  )

patient_feature_df <- lm_patient_df %>%
  full_join(mono_patient_df, by = "patient_id") %>%
  full_join(lymphoid_score_patient_df, by = "patient_id") %>%
  full_join(mono_score_patient_df, by = "patient_id") %>%
  arrange(patient_id)

write.table(
  patient_feature_df,
  file = "data/20260309_pilot/results/version1/data_export/20260320_patient_level_immune_features.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

# 11. raw_query metadata로 cell_meta 생성 + join
raw_query <- readRDS('data/20260309_pilot/results/20260313_raw_query_annotation_attached.rds')
#DefaultAssay(raw_query) <- "RNA"
#raw_query <- NormalizeData(raw_query)

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
  select(
    cell_id,
    patient_id,
    source,
    collection_batch,
    timepoint,
    celltype_main,
    nCount_RNA,
    Sample.ID,
    Study.ID,
    Age_IO,
    sex,
    Histology,
    smoking,
    ecog,
    drug,
    PD_event,
    PFS,
    Death_Event,
    OS,
    EGFR_Mutation_Status,
    IO_Line,
    Previous_palliative_chemo,
    Previous_palliative_target,
    PD.L1_TPS
  ) %>%
  left_join(cell_score_df, by = "cell_id") %>%
  left_join(patient_feature_df, by = "patient_id")

arrow::write_parquet(
  cell_meta,
  "data/20260309_pilot/results/version1/data_export/20260320_cell_meta.parquet"
)

###
DimPlot(seu_integrated,reduction = "umap",  group.by='manual.cluster',label=T) #20260312_integrated_pca20_manualcluster_umap
# hard to distinguish cytotoxic and activation /naive and helper-memory like
# hard to find exhaustive cluster 
FeaturePlot(seu_integrated, features = c('NKG7', 'CCL5', 'KLRD1', 'KLRF1','FGFBP2', 'PRF1', 'CTSW', 'GNLY','GZMB', 'GZMH','CCL4')) #cytotoxic #8x12 #20260320_cytotoxic_gene_set
FeaturePlot(seu_integrated, features = c('CD69','IL2RA','TNFRSF9','IFNG','IL2','TNF')) #activated # no enrichemnt
FeaturePlot(seu_integrated, features = c('TCF7', 'LEF1', 'CCR7', 'IL7R','MAL', 'LTB')) #naive/stem-like #6x6 #20260320_naive_gene_set
FeaturePlot(seu_integrated, features = c('IL7R', 'LTB','MAL', 'IL6R')) #helper-like/memory-like #4x6  #20260320_helper_gene_set
FeaturePlot(seu_integrated, features = c('PDCD1', 'HAVCR2', 'TIGIT', 'LAG3','CTLA4', 'TOX','CXCL13','ENTPD1','LAYN')) #exhaustion <- no enrichment
navie = c('TCF7', 'LEF1', 'CCR7', 'IL7R','MAL', 'LTB')
cytotoxic = c('NKG7', 'CCL5', 'KLRD1', 'KLRF1','FGFBP2', 'PRF1', 'CTSW', 'GNLY','GZMB', 'GZMH','CCL4')

# classic - non-classic
FeaturePlot(seu_integrated, features = c('CD14', 'VCAN', 'S100A8', 'S100A9', 'LYZ')) #CD14, classic #6x6 #20260320_CD14_mono_gene_set
FeaturePlot(seu_integrated, features = c('FCGR3A', 'LST1', 'CDKN1C', 'MS4A7', 'IFITM3', 'CX3CR1','GPBAR1','LRRC25', 'CALHM6','CSF1R', 'MS4A7', 'LILRB2', 'MAFB')) #CD14, classic #8x12 #20260320_CD16_mono_gene_set

cell_df <- seu_integrated@meta.data[seu_integrated@meta.data$time =='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  dplyr::select(sample, manual.cluster)
count_df <- seu_integrated@meta.data[seu_integrated@meta.data$time =='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  count(sample, manual.cluster, name = "n_cells")
total_df <- seu_integrated@meta.data[seu_integrated@meta.data$time =='base' & seu_integrated@meta.data$cohort =='Lee_p1',] %>%
  count(sample, name = "total_cells")
prop_df <- count_df %>%
  left_join(total_df, by = "sample") %>%
  mutate(prop = n_cells / total_cells)




