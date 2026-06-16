library(data.table)
library(ggplot2)
library(reshape2)
library(kernlab)
library(Rtsne)
library(uwot)
library(FNN)
library(igraph)
library(readxl)
install.packages(c("FactoMineR", "factoextra"))
library(FactoMineR)
library(factoextra)

### 
# 0. meta
###
meta <- fread('data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt')
meta$sample_id <- gsub('_','-',meta$sample_id)
dim(meta)
##########
## 1. load full feature data
###########
B <- fread('data/20260309_pilot/results/version2/B/final_feature_df_with_cellrank_b3_filt.csv')
CD4T <- fread('data/20260309_pilot/results/version2/CD4T/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv')
CD8T <- fread('data/20260309_pilot/results/version2/CD8T/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv')
Monocyte <- fread('data/20260309_pilot/results/version2/Monocyte/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv')
NK <- fread('data/20260309_pilot/results/version2/NK/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv')
nonconventional_T <- fread('data/20260309_pilot/results/version2/nonconventional_T/CellRank_output_from_clean_subset/final_feature_df_with_cellrank_b3_filt.csv')

lst <- list(
  B = B,
  CD4T = CD4T,
  CD8T = CD8T,
  Monocyte = Monocyte,
  NK = NK,
  nonconventional_T = nonconventional_T
)

lst2 <- Map(function(dt, nm) {
  dt <- copy(dt)
  cols_to_rename <- setdiff(names(dt), "sample")
  setnames(dt, cols_to_rename, paste0(nm, "__", cols_to_rename))
  dt
}, lst, names(lst))

merged_dt <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE), lst2)

merged_dt <- merged_dt[merged_dt$sample %in% meta$sample_id,]
##########
## A-2. pick features included for patiaent clustering
###########
ov_rank <- fread('data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top10.csv')
ov_rank <- fread('data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top20.csv')
ov_rank <- fread('data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer/outer_search_validation/repeated_top_features__top30.csv')
head(ov_rank[c(order(ov_rank$mean_group_rank)),c('representative_feature_name','n_folds_in_top_n','mean_raw_rank','mean_group_rank','group_rank_fold_01','group_rank_fold_02','group_rank_fold_03','roc_auc_mean_fold_01','roc_auc_mean_fold_02','roc_auc_mean_fold_03')])
ov_rank_2f <- ov_rank[ov_rank$n_folds_in_top_n >=2,]
dim(ov_rank_2f)
ov_rank <- ov_rank[c(order(ov_rank$mean_group_rank)),]
ggplot(ov_rank, aes(x=factor(representative_feature_name,levels=ov_rank$representative_feature_name), y=mean_group_rank))+geom_point()
ov_rank_melt <- melt(ov_rank, id.vars = c('representative_feature_name','n_folds_in_top_n'), measure.vars = c('mean_group_rank','group_rank_fold_01','group_rank_fold_02','group_rank_fold_03'))
ggplot(ov_rank_melt[ov_rank_melt$variable !='mean_group_rank',], aes(x=factor(representative_feature_name,levels=ov_rank$representative_feature_name), y=value, col=variable,group=variable))+geom_line()+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=2))#5x10 #20260427_version2_trial1_fold_group_rank

#######
# A-3. leverage results from endpoint prediction
#######
top30 <- ov_rank[seq(1,30),]$representative_feature_name
cols_keep <- c('sample',top30)
merged_dt_filt <- merged_dt[,..cols_keep]
merged_dt_filt.melt <- melt(merged_dt_filt, id.vars='sample',measure.vars=top30)
ggplot(merged_dt_filt.melt,aes(x=value, col=variable))+geom_density()+facet_wrap(.~variable, scales='free')+guides(col='none')

ov_rank_tmp <- ov_rank[ov_rank$representative_feature_name %in% top30,]
ov_rank_tmp$baseline <- ov_rank_tmp$mean_roc_auc_mean - ov_rank_tmp$mean_delta_roc_auc_mean
ov_rank_tmp$baseline_fold1 <- ov_rank_tmp$roc_auc_mean_fold_01 - ov_rank_tmp$delta_roc_auc_mean_fold_01
ov_rank_tmp$baseline_fold2 <- ov_rank_tmp$roc_auc_mean_fold_02 - ov_rank_tmp$delta_roc_auc_mean_fold_02
ov_rank_tmp$baseline_fold3 <- ov_rank_tmp$roc_auc_mean_fold_03 - ov_rank_tmp$delta_roc_auc_mean_fold_03
ov_rank_tmp <- ov_rank_tmp[c(order(ov_rank_tmp$mean_roc_auc_mean)),]
ov_rank_tmp_melt <- melt(ov_rank_tmp, id.vars=c('representative_feature_name'), measure.vars = c('baseline_fold1','baseline_fold2','baseline_fold3','roc_auc_mean_fold_01','roc_auc_mean_fold_02','roc_auc_mean_fold_03'))
ov_rank_tmp_melt <- as.data.table(ov_rank_tmp_melt)
ov_rank_tmp_melt[, fold := fifelse(
  grepl("fold1|fold_01", variable), "fold1",
  fifelse(
    grepl("fold2|fold_02", variable), "fold2",
    fifelse(grepl("fold3|fold_03", variable), "fold3", NA_character_)
  )
)]

ov_rank_tmp_melt[, line_type := fifelse(
  variable %in% c("baseline_fold1", "baseline_fold2", "baseline_fold3"),
  "baseline",
  "model"
)]
ggplot(
  ov_rank_tmp_melt,
  aes(
    x = factor(
      representative_feature_name,
      levels = ov_rank_tmp$representative_feature_name
    ),
    y = value,
    color = fold,
    group = variable,
    linetype = line_type
  )
) +
  geom_line() +
  geom_point() +
  scale_color_manual(
    values = c(
      fold1 = "grey10",
      fold2 = "grey50",
      fold3 = "grey80"
    )
  ) +
  scale_linetype_manual(
    values = c(
      baseline = "dashed",
      model = "solid"
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Feature",
    y = "ROC-AUC",
    color = "Fold",
    linetype = ""
  )#8x12 #20260429_outervalidation_fold_rocauc
#####
# A-5.결측값 대처 (median imputation)
#####
colSums(is.na(merged_dt_filt))
for (col in top30) {
  merged_dt_filt[[col]] <- as.numeric(merged_dt_filt[[col]])
  
  if (anyNA(merged_dt_filt[[col]])) {
    med_val <- median(merged_dt_filt[[col]], na.rm = TRUE)
    merged_dt_filt[[col]][is.na(merged_dt_filt[[col]])] <- med_val
  }
}
colSums(is.na(merged_dt_filt))
#####
# A-6.robust min-max normalization
#####
robust_minmax <- function(x, q_low = 0.02, q_high = 0.98) {
  x <- as.numeric(x)
  
  qs <- quantile(
    x,
    probs = c(q_low, q_high),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  
  q02 <- qs[1]
  q98 <- qs[2]
  
  if (is.na(q02) || is.na(q98) || q98 <= q02) {
    return(rep(0, length(x)))
  }
  
  z <- (x - q02) / (q98 - q02)
  z <- pmin(pmax(z, 0), 1)
  return(z)
}
df_norm <- copy(merged_dt_filt)

for (col in top30) {
  df_norm[[col]] <- robust_minmax(df_norm[[col]])
}

df_norm.melt <- melt(df_norm, id.vars='sample',measure.vars=top30)
ggplot(df_norm.melt,aes(x=value, col=variable))+geom_density()+facet_wrap(.~variable, scales='free')+guides(col='none')

summary(df_norm[, ..top30])
######
# A-7. clustering용 matrix 만들기, k clustering
######
x_mat <- as.matrix(df_norm[, ..top30])
rownames(x_mat) <- df_norm[['sample']]

####
# A-FAMD version
# move to line 1745
#####
set.seed(123)

k_range <- 2:8

k_summary_list <- list()
k_pairwise_list <- list()
k_cluster_list <- list()

for (k in k_range) {
  set.seed(123 + k)
  
  sc_res <- specc(
    x = x_mat,
    centers = k
  )
  
  clusters <- as.integer(sc_res)
  cluster_pairs <- combn(sort(unique(clusters)), 2)
  
  pairwise_res_list <- list()
  
  for (i in seq_len(ncol(cluster_pairs))) {
    a <- cluster_pairs[1, i]
    b <- cluster_pairs[2, i]
    
    idx_a <- which(clusters == a)
    idx_b <- which(clusters == b)
    
    pvals <- sapply(top30, function(feat) {
      tryCatch(
        wilcox.test(
          x_mat[idx_a, feat],
          x_mat[idx_b, feat],
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      )
    })
    
    padj <- p.adjust(pvals, method = "bonferroni")
    n_sig <- sum(padj < 0.05, na.rm = TRUE)
    
    pairwise_res_list[[i]] <- data.table(
      k = k,
      cluster_a = a,
      cluster_b = b,
      n_sig_features = n_sig
    )
  }
  
  pairwise_res <- rbindlist(pairwise_res_list)
  
  k_summary_list[[as.character(k)]] <- data.table(
    k = k,
    median_sig_features = median(pairwise_res$n_sig_features, na.rm = TRUE),
    mean_sig_features = mean(pairwise_res$n_sig_features, na.rm = TRUE),
    min_sig_features = min(pairwise_res$n_sig_features, na.rm = TRUE),
    max_sig_features = max(pairwise_res$n_sig_features, na.rm = TRUE),
    min_cluster_size = min(table(clusters)),
    max_cluster_size = max(table(clusters))
  )
  
  k_pairwise_list[[as.character(k)]] <- pairwise_res
  
  k_cluster_list[[as.character(k)]] <- data.table(
    sample = rownames(x_mat),
    k = k,
    cluster = clusters
  )
}

k_summary_dt <- rbindlist(k_summary_list)
k_summary_dt[order(k)]

k_summary_dt <- k_summary_dt[order(k)]

best_score <- max(k_summary_dt$median_sig_features, na.rm = TRUE)
best_k_candidates <- k_summary_dt[median_sig_features == best_score, k]

best_score
best_k_candidates

ggplot(k_summary_dt, aes(x = k, y = median_sig_features)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    x = "Number of clusters (K)",
    y = "Median number of significant features",
    title = "K selection based on pairwise feature separation"
  )#3x4 #20260427_version2_trial1_k
best_k <- max(best_k_candidates)
best_k =2
second_best_k =4
# final_cluster_dt <- k_cluster_list[[as.character(best_k)]]
# final_pairwise_dt <- k_pairwise_list[[as.character(best_k)]]
final_cluster_dt <- k_cluster_list[[as.character(second_best_k)]]
final_pairwise_dt <- k_pairwise_list[[as.character(second_best_k)]]

final_cluster_dt[1:10]
final_pairwise_dt

######
# A-8. visualization
#####
plot_dt <- copy(final_cluster_dt)
plot_dt[, cluster := factor(cluster)]
head(plot_dt)
plot_dt <- merge(plot_dt, meta, by.x='sample', by.y='sample_id')
# pca
pca_res <- prcomp(x_mat, center = TRUE, scale. = FALSE)

var_explained <- pca_res$sdev^2
pct_var <- var_explained / sum(var_explained) * 100
cum_var <- cumsum(pct_var)

pca_var_dt <- data.table(
  PC = seq_along(pct_var),
  pct_var = pct_var,
  cum_var = cum_var
)

ggplot(pca_var_dt, aes(x = PC, y = pct_var)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  labs(
    x = "Principal component (PC)",
    y = "Percent variance explained",
    title = "PCA elbow plot"
  )#3x4 #20260427_version2_trial1_pc_variance_explained


pca_dt <- data.table(
  sample = rownames(x_mat),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3],
  PC4 = pca_res$x[, 4],
  PC5 = pca_res$x[, 5],
  PC6 = pca_res$x[, 6],
  PC7 = pca_res$x[, 7],
  PC8 = pca_res$x[, 8],
  PC9 = pca_res$x[, 9],
  PC10 = pca_res$x[, 10],
  PC11 = pca_res$x[, 11],
  PC12 = pca_res$x[, 12]
)

pca_plot_dt <- merge(pca_dt, plot_dt, by = "sample")
pca_plot_dt_batch <- merge(pca_plot_dt, old[,c('patient_id','batch')],by.x='sample',by.y='patient_id')
# ggplot(pca_plot_dt, aes(x = PC2, y = PC4, color = cluster)) +
#   geom_point(size = 3) +
#   theme_bw()
# table(pca_plot_dt$cluster,pca_plot_dt$`Binarized response`)
# 
ggplot(pca_plot_dt_batch, aes(x = PC1, y = PC2, color = grepl('^M',batch), label = batch,shape=cluster)) +
  geom_point(size = 3) +
  geom_text_repel(vjust = -0.5, size = 3) +
  theme_bw()#20260427_version2_trial1_pca12_cluster #20260427_version2_trial1_pca23_cluster #20260427_version2_trial1_pca12_cluster_v2
# table(pca_plot_dt$cluster,pca_plot_dt$`Binarized response`)

ggplot(pca_plot_dt_batch, aes(x = PC2, y = PC3, color = cluster, label = batch, shape=cluster)) +
  geom_point(size =3) +
  geom_text(vjust = -1, size = 2) +
  theme_bw()#20260429_version2_trial1_pc23_cluster
ggplot(pca_plot_dt_batch, aes(x = PC2, y = PC3, color = `Binarized response`, label = batch, shape=cluster)) +
  geom_point(size =3) +
  geom_text(vjust = -1, size = 2) +
  theme_bw()#5x6 #20260427_version2_trial1_pca12 #20260427_version2_trial1_pca23 #20260427_version2_trial1_pca24 #20260429_version2_trial1_pc23_binarizedresponse

pca_plot_dt_batch$PFS_180 <- pmin(pca_plot_dt_batch$`PFS (Days)`, 180)
ggplot(pca_plot_dt_batch, aes(x = PC2, y = PC4, color = PFS_180, label = batch, shape=cluster)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 2) +
  theme_bw()#5x6 #20260427_version2_trial1_pca12 #20260427_version2_trial1_pca23 #20260427_version2_trial1_pca24 #20260429_version2_trial1_pc23_PFS

ggplot(pca_plot_dt_batch, aes(x = PC2, y = PC4, color = `RECIST response`, label = batch, shape=cluster)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 2) +
  theme_bw()#20260429_version2_trial1_pc23_RECIST

### 어떤 PC 조합이 어떤 metadata 정보를 잘 설명하는가
## 1.전체 PC1–PC8이 각 metadata를 얼마나 설명하는지
library(vegan)
dt <- copy(pca_plot_dt_batch)

pc_cols <- paste0("PC", 1:9)

# ID나 이미 PCA에서 나온 column은 제외
exclude_cols <- c(
  "sample", "Study ID", "Sample ID",
  pc_cols
)

meta_cols <- setdiff(names(dt), exclude_cols)

# 값이 2개 이상 있는 column만 사용
meta_cols <- meta_cols[
  sapply(dt[, ..meta_cols], function(x) data.table::uniqueN(na.omit(x)) >= 2)
]

run_adonis_one <- function(dt, var, pcs = pc_cols, nperm = 999) {
  sub <- dt[, c(pcs, var), with = FALSE]
  sub <- na.omit(sub)
  
  if (nrow(sub) < 5) return(NULL)
  if (data.table::uniqueN(sub[[var]]) < 2) return(NULL)
  
  # PC는 numeric 보장
  for (p in pcs) sub[[p]] <- as.numeric(sub[[p]])
  
  # metadata column을 안전한 이름으로 복사
  if (is.character(sub[[var]]) || is.factor(sub[[var]])) {
    sub[, meta_var := factor(get(var))]
  } else {
    sub[, meta_var := as.numeric(get(var))]
  }
  
  X <- as.matrix(sub[, ..pcs])
  d <- dist(X)
  
  fit <- vegan::adonis2(
    d ~ meta_var,
    data = as.data.frame(sub),
    permutations = nperm
  )
  
  data.table(
    var = var,
    pc_set = paste(pcs, collapse = "+"),
    n = nrow(sub),
    n_level = data.table::uniqueN(sub$meta_var),
    R2 = fit$R2[1],
    F = fit$F[1],
    p = fit$`Pr(>F)`[1]
  )
}

global_pc_res <- rbindlist(
  lapply(meta_cols, function(v) run_adonis_one(dt, v, pcs = pc_cols)),
  fill = TRUE
)

global_pc_res[, p_adj := p.adjust(p, method = "BH")]
setorder(global_pc_res, -R2)

global_pc_res
global_pc_res[order(-R2)]
## 2.어떤 PC 하나가 어떤 metadata를 잘 설명하는지 heatmap
single_pc_res <- rbindlist(
  lapply(meta_cols, function(v) {
    rbindlist(
      lapply(pc_cols, function(pc) {
        run_adonis_one(dt, v, pcs = pc)
      }),
      fill = TRUE
    )
  }),
  fill = TRUE
)

single_pc_res[, PC := pc_set]
single_pc_res[, p_adj := p.adjust(p, method = "BH")]
single_pc_res[, neglog10_padj := -log10(p_adj)]

single_pc_res[order(-R2)]
ggplot(single_pc_res, aes(x = PC, y = reorder(var, R2), fill = R2)) +
  geom_tile() +
  theme_bw() +
  labs(
    x = "PC",
    y = "Metadata column",
    fill = "R2",
    title = "How strongly each PC is associated with each metadata column"
  )#5x6 #20260427_version2_trial1 #20260427_version2_trial1_exceptBatchCluster

# use line 1298 function
res_pc3 <- plot_pc_feature_heatmap(
  pca_res = pca_res,
  x_input = x_mat,
  pc = "PC3",
  top_n = 5,
  meta_dt = pca_plot_dt_batch,
  sample_id_col = "sample",
  annotation_cols = c("cluster", "ECOG _PS"),
  scale_by_feature = FALSE
)#3x13 #20260427_version2_trial1_PC2 #20260427_version2_trial1_PC3

meta_check <- unique(
  as.data.table(pca_plot_dt_batch)[, .(
    sample,
    cluster,
    `Binarized response`,
    `RECIST response`
  )],
  by = "sample"
)

colSums(is.na(meta_check))

#####
# tsne
set.seed(123)

tsne_res <- Rtsne(
  x_mat,
  dims = 2,
  perplexity = min(10, floor((nrow(x_mat) - 1) / 3)),
  check_duplicates = FALSE,
  pca = FALSE,
  max_iter = 1000,
  verbose = TRUE
)

tsne_dt <- data.table(
  sample = rownames(x_mat),
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2]
)

tsne_plot_dt <- merge(tsne_dt, plot_dt, by = "sample")

ggplot(tsne_plot_dt, aes(x = tSNE1, y = tSNE2, color = `Binarized response`, label = sample)) +
  geom_point(size = 3) +geom_text(vjust = -0.5, size = 3) +
  theme_bw()

ggplot(tsne_plot_dt, aes(x = tSNE1, y = tSNE2, color = `Drug`, label = sample)) +
  geom_point(size = 3) +
  theme_bw()

tsne_plot_dt.melt <- melt(tsne_plot_dt, id.vars=c('sample','tSNE1','tSNE2','k'),measure.vars=c('cluster','Binarized response','RECIST response','Sex M: 1  F: 2','Age_IO','Histology_mod','Smoking Never: 0 Ex: 1 Current: 2','ECOG _PS','Drug','PD_Event','PFS (Days)','OS (Days)','EGFR_Mutation_Status_mod','IO_Line','Previous_palliative_chemo','Previous_palliative_target','PD-L1_TPS_mod'))

ggplot(tsne_plot_dt.melt, aes(x = tSNE1, y = tSNE2, color = value, label = sample)) +
  geom_point(size = 1) +facet_wrap(.~variable)+guides(col='none')+
  theme_bw()

# umap
set.seed(123)
#pc_use <- c(4,5,6,7,8,9,10)

umap_input <- pca_res$x[, 2:12, drop = FALSE]

umap_mat <- uwot::umap(
  umap_input,
  n_neighbors = 10,
  min_dist = 0.3,
  metric = "euclidean"
)

umap_dt <- data.table(
  sample = rownames(pca_res$x),
  UMAP1 = umap_mat[, 1],
  UMAP2 = umap_mat[, 2]
)

umap_plot_dt <- merge(
  umap_dt,
  pca_plot_dt_batch,
  by = "sample"
)

ggplot(
  umap_plot_dt,
  aes(x = UMAP1, y = UMAP2, color = cluster,label = batch)
) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw()#20260427_version2_trial1_umap

ggplot(
  umap_plot_dt,
  aes(x = UMAP1, y = UMAP2, color = `Binarized response`,label = batch)
) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw()#20260427_version2_trial1_umap_v2

# umap_res <- umap(
#   x_mat,
#   n_neighbors = min(10, nrow(x_mat) - 1),
#   min_dist = 0.3,
#   metric = "euclidean"
# )
# 
# umap_dt <- data.table(
#   sample = rownames(x_mat),
#   UMAP1 = umap_res[, 1],
#   UMAP2 = umap_res[, 2]
# )

# umap_plot_dt <- merge(umap_dt, plot_dt, by = "sample", all.x = TRUE)




#######
# B. use higly variable features (not using endpoint)
# B-2. make feature matrix
######
merged_dt_anno <- merge(merged_dt, meta, by.x='sample', by.y='sample_id')
non_feature_cols <-colnames(meta)
non_feature_cols <- non_feature_cols[non_feature_cols !='sample_id']
feature_cols <- colnames(merged_dt)
feature_cols <- feature_cols[feature_cols !='sample']

raw_df <- copy(merged_dt_anno)
raw_df[, row_id := .I]

sample_id <- raw_df$sample
row_id <- raw_df$row_id

x_df <- raw_df[, ..feature_cols]

# df_use <- copy(merged_dt_anno)
# # sample ID
# sample_id <- df_use$sample
# numeric feature matrix
# x_df <- df_use[, ..feature_cols]

#####
# B-3. 결측값 처리 (median imputation)
######
colSums(is.na(x_df))[colSums(is.na(x_df)) > 0]
for (j in seq_len(ncol(x_df))) {
  if (anyNA(x_df[[j]])) {
    x_df[[j]][is.na(x_df[[j]])] <- median(x_df[[j]], na.rm = TRUE)
  }
}
colSums(is.na(x_df))[colSums(is.na(x_df)) > 0]

####
# B-4. remove near zero variance
####
feature_var <- sapply(x_df, var, na.rm = TRUE)
feature_sd  <- sapply(x_df, sd, na.rm = TRUE)
feature_nuniq <- sapply(x_df, function(v) length(unique(v)))

nzv_dt <- data.table(
  feature = names(x_df),
  variance = feature_var,
  sd = feature_sd,
  n_unique = feature_nuniq
)

head(nzv_dt[order(variance)])

keep_nzv <- nzv_dt[
  variance > 1e-8 & n_unique > 2,
  feature
]
length(keep_nzv) #944
x_nzv <- x_df[, ..keep_nzv]

###
# B-4.2. feature pruning
# B-4.2.1 helper: feature name parsing
###
known_prefixes <- c(
  "Centroid_", "curated_gene_", "denovo_gene_",
  "ilr_", "cellrank_", "cellrank_cr_b3"
)

starts_with_any <- function(x, prefixes) {
  any(startsWith(x, prefixes))
}

split_feature_parts <- function(feature_name) {
  feature <- as.character(feature_name)
  stem <- ""
  core <- feature
  
  if (grepl("__", feature, fixed = TRUE)) {
    sp <- strsplit(feature, "__", fixed = TRUE)[[1]]
    left <- sp[1]
    right <- paste(sp[-1], collapse = "__")
    
    if (!starts_with_any(left, known_prefixes)) {
      stem <- left
      core <- right
    }
  }
  
  list(feature = feature, stem = stem, core = core)
}

infer_feature_family <- function(feature_name) {
  parts <- split_feature_parts(feature_name)
  feature <- parts$feature
  stem <- parts$stem
  core <- parts$core
  
  if (startsWith(core, "Centroid_")) {
    return(if (nzchar(stem)) paste0("centroid::", stem) else "centroid")
  }
  
  if (startsWith(core, "curated_gene_")) {
    suffixes <- c("__eigengene", "__pc1", "__pc2", "__singscore", "__ssgsea2_es")
    signature_stem <- core
    for (suffix in suffixes) {
      if (endsWith(core, suffix)) {
        signature_stem <- substr(core, 1, nchar(core) - nchar(suffix))
        break
      }
    }
    return(
      if (nzchar(stem)) {
        paste0("signature::", stem, "::", signature_stem)
      } else {
        paste0("signature::", signature_stem)
      }
    )
  }
  
  if (startsWith(core, "cellrank_")) {
    return(if (nzchar(stem)) paste0("cellrank::", stem) else "cellrank")
  }
  
  if (nzchar(stem)) {
    return(paste0("file::", stem))
  }
  
  if (startsWith(feature, "PC")) {
    return("pc")
  }
  
  "misc"
}

infer_family_cap_group <- function(feature_name) {
  parts <- split_feature_parts(feature_name)
  stem <- parts$stem
  core <- parts$core
  
  if (!startsWith(core, "curated_gene_")) {
    return(NA_character_)
  }
  
  suffixes <- c("__eigengene", "__pc1", "__pc2", "__singscore", "__ssgsea2_es")
  signature_stem <- core
  for (suffix in suffixes) {
    if (endsWith(core, suffix)) {
      signature_stem <- substr(core, 1, nchar(core) - nchar(suffix))
      break
    }
  }
  
  if (nzchar(stem)) {
    paste0("signature::", stem, "::", signature_stem)
  } else {
    paste0("signature::", signature_stem)
  }
}

infer_feature_level <- function(feature_name) {
  feature <- as.character(feature_name)
  core <- feature
  
  if (grepl("__", feature, fixed = TRUE)) {
    sp <- strsplit(feature, "__", fixed = TRUE)[[1]]
    left <- sp[1]
    right <- paste(sp[-1], collapse = "__")
    core <- if (starts_with_any(left, known_prefixes)) feature else right
  }
  
  if (startsWith(core, "Centroid_")) return("PCA_centroid_continuum")
  if (startsWith(core, "curated_gene_") && endsWith(core, "__ssgsea2_es")) return("pseudobulk_curated_ssgsea2")
  if (startsWith(core, "curated_gene_") && endsWith(core, "__singscore")) return("pseudobulk_curated_singscore")
  if (startsWith(core, "curated_gene_") && endsWith(core, "__pc1")) return("pseudobulk_curated_pc1")
  if (startsWith(core, "curated_gene_") && endsWith(core, "__pc2")) return("pseudobulk_curated_pc2")
  if (startsWith(core, "curated_gene_") && endsWith(core, "__eigengene")) return("pseudobulk_curated_eigengene")
  if (startsWith(core, "denovo_gene_")) return("pseudobulk_de_novo")
  if (startsWith(core, "ilr_")) return("ILR_composition")
  if (startsWith(core, "cellrank_cr_b3")) return("status_dynamics")
  
  "patient_feature"
}

###
# B-4.2.2. helper: variablity, normality
###

qq_normality_score <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) < 3) return(NA_real_)
  if (length(unique(x)) < 3) return(-Inf)
  
  x <- sort(x)
  theo <- qnorm(ppoints(length(x)))
  
  out <- suppressWarnings(cor(x, theo, method = "pearson"))
  as.numeric(out)
}

calc_feature_metrics <- function(dt) {
  dt <- as.data.table(dt)
  
  res <- lapply(names(dt), function(feat) {
    x <- as.numeric(dt[[feat]])
    
    data.table(
      feature = feat,
      variance = var(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      mad = mad(x, na.rm = TRUE),
      iqr = IQR(x, na.rm = TRUE),
      n_unique = length(unique(x)),
      normality_score = qq_normality_score(x)
    )
  })
  
  rbindlist(res)
}
###
# B-4.2.3.helper: representative feature 고르기
###
choose_representative_feature <- function(metrics_sub,
                                          variability_metric = c("mad", "variance", "sd", "iqr"),
                                          preference = c("normality_then_variability", "variability_then_normality")) {
  variability_metric <- match.arg(variability_metric)
  preference <- match.arg(preference)
  
  metrics_sub <- as.data.table(copy(metrics_sub))
  metrics_sub[, variability_value := get(variability_metric)]
  
  if (preference == "normality_then_variability") {
    setorder(metrics_sub, -normality_score, -variability_value, -n_unique, feature)
  } else {
    setorder(metrics_sub, -variability_value, -normality_score, -n_unique, feature)
  }
  
  metrics_sub$feature[1]
}

choose_better_of_two <- function(f1, f2, metrics_dt,
                                 variability_metric = "mad",
                                 preference = "normality_then_variability") {
  sub <- metrics_dt[feature %in% c(f1, f2)]
  choose_representative_feature(sub,
                                variability_metric = variability_metric,
                                preference = preference)
}

###
# B-4.2.4 make feature metadata table
###
make_feature_metadata <- function(feature_names) {
  data.table(
    feature = feature_names,
    family = vapply(feature_names, infer_feature_family, character(1)),
    family_cap_group = vapply(feature_names, infer_family_cap_group, character(1)),
    feature_level = vapply(feature_names, infer_feature_level, character(1))
  )
}

###
# B-4.2.5 family cap pruning
###
family_cap_prune <- function(dt,
                             metrics_dt,
                             feature_meta,
                             group_col = "family_cap_group",
                             variability_metric = "mad",
                             preference = "normality_then_variability") {
  dt <- as.data.table(dt)
  metrics_dt <- as.data.table(metrics_dt)
  feature_meta <- as.data.table(feature_meta)
  
  meta <- merge(feature_meta, metrics_dt, by = "feature", all.x = TRUE)
  
  grouped <- meta[!is.na(get(group_col)) & get(group_col) != ""]
  ungrouped <- meta[is.na(get(group_col)) | get(group_col) == ""]
  
  keep_features <- c()
  drop_log <- list()
  
  # grouped features: representative 1개 선택
  if (nrow(grouped) > 0) {
    for (grp in unique(grouped[[group_col]])) {
      sub <- grouped[get(group_col) == grp]
      
      if (nrow(sub) == 1) {
        keep_features <- c(keep_features, sub$feature)
      } else {
        rep_feat <- choose_representative_feature(
          sub,
          variability_metric = variability_metric,
          preference = preference
        )
        keep_features <- c(keep_features, rep_feat)
        
        dropped <- setdiff(sub$feature, rep_feat)
        if (length(dropped) > 0) {
          drop_log[[length(drop_log) + 1]] <- data.table(
            step = "family_cap",
            group_name = grp,
            kept_feature = rep_feat,
            dropped_feature = dropped
          )
        }
      }
    }
  }
  
  # ungrouped features: 그대로 유지
  keep_features <- c(keep_features, ungrouped$feature)
  keep_features <- unique(keep_features)
  
  drop_dt <- if (length(drop_log) > 0) rbindlist(drop_log) else data.table()
  
  list(
    kept_features = keep_features,
    kept_dt = dt[, ..keep_features],
    drop_log = drop_dt
  )
}

###
# B-4.2.6 pairwise correlation pruning
###
correlation_prune <- function(dt,
                              metrics_dt,
                              threshold = 0.8,
                              method = c("pearson", "spearman"),
                              variability_metric = "mad",
                              preference = "normality_then_variability") {
  method <- match.arg(method)
  dt <- as.data.table(dt)
  metrics_dt <- as.data.table(metrics_dt)
  
  remaining <- names(dt)
  
  cor_mat <- suppressWarnings(
    cor(as.matrix(dt), use = "pairwise.complete.obs", method = method)
  )
  
  drop_log <- list()
  
  repeat {
    current_cor <- cor_mat[remaining, remaining, drop = FALSE]
    diag(current_cor) <- NA_real_
    
    max_abs_cor <- suppressWarnings(max(abs(current_cor), na.rm = TRUE))
    
    if (!is.finite(max_abs_cor) || max_abs_cor < threshold) {
      break
    }
    
    hit <- which(abs(current_cor) == max_abs_cor, arr.ind = TRUE)[1, ]
    f1 <- rownames(current_cor)[hit[1]]
    f2 <- colnames(current_cor)[hit[2]]
    
    keep_feat <- choose_better_of_two(
      f1, f2, metrics_dt,
      variability_metric = variability_metric,
      preference = preference
    )
    drop_feat <- setdiff(c(f1, f2), keep_feat)[1]
    
    drop_log[[length(drop_log) + 1]] <- data.table(
      step = "correlation_prune",
      feature_a = f1,
      feature_b = f2,
      abs_correlation = abs(current_cor[f1, f2]),
      kept_feature = keep_feat,
      dropped_feature = drop_feat
    )
    
    remaining <- setdiff(remaining, drop_feat)
  }
  
  drop_dt <- if (length(drop_log) > 0) rbindlist(drop_log) else data.table()
  
  list(
    kept_features = remaining,
    kept_dt = dt[, ..remaining],
    drop_log = drop_dt
  )
}

###
# B-4.2.7 VIF pruning
###

calc_vif_table <- function(dt) {
  dt <- as.data.table(dt)
  x <- as.matrix(dt)
  
  out <- lapply(seq_len(ncol(x)), function(j) {
    y <- x[, j]
    X <- x[, -j, drop = FALSE]
    
    if (ncol(X) == 0) {
      return(data.table(feature = colnames(x)[j], vif = 1))
    }
    
    fit <- lm.fit(x = cbind(Intercept = 1, X), y = y)
    
    rss <- sum(fit$residuals^2)
    tss <- sum((y - mean(y))^2)
    
    if (!is.finite(tss) || tss <= .Machine$double.eps) {
      vif_val <- Inf
    } else {
      r2 <- 1 - rss / tss
      r2 <- min(max(r2, 0), 0.999999999999)
      vif_val <- 1 / (1 - r2)
    }
    
    data.table(feature = colnames(x)[j], vif = vif_val)
  })
  
  rbindlist(out)
}

vif_prune <- function(dt,
                      metrics_dt,
                      vif_threshold = 10,
                      variability_metric = "mad",
                      preference = "normality_then_variability",
                      verbose = TRUE) {
  dt <- as.data.table(dt)
  metrics_dt <- as.data.table(metrics_dt)
  
  remaining <- names(dt)
  drop_log <- list()
  iter <- 0
  
  repeat {
    iter <- iter + 1
    
    if (length(remaining) <= 1) break
    
    vif_dt <- calc_vif_table(dt[, ..remaining])
    max_vif <- max(vif_dt$vif, na.rm = TRUE)
    
    if (verbose) {
      message("VIF iteration ", iter, ": max VIF = ", round(max_vif, 4),
              " | n_features = ", length(remaining))
    }
    
    if (!is.finite(max_vif) || max_vif > vif_threshold) {
      candidates <- vif_dt[vif == max_vif, feature]
      
      if (length(candidates) == 1) {
        drop_feat <- candidates[1]
      } else {
        # tie면 "덜 좋은" feature를 drop
        sub <- metrics_dt[feature %in% candidates]
        keep_best <- choose_representative_feature(
          sub,
          variability_metric = variability_metric,
          preference = preference
        )
        drop_feat <- setdiff(candidates, keep_best)[1]
      }
      
      drop_log[[length(drop_log) + 1]] <- data.table(
        step = "vif_prune",
        iteration = iter,
        max_vif = max_vif,
        dropped_feature = drop_feat
      )
      
      remaining <- setdiff(remaining, drop_feat)
    } else {
      break
    }
  }
  
  final_vif <- if (length(remaining) >= 2) {
    calc_vif_table(dt[, ..remaining])
  } else {
    data.table(feature = remaining, vif = 1)
  }
  
  drop_dt <- if (length(drop_log) > 0) rbindlist(drop_log) else data.table()
  
  list(
    kept_features = remaining,
    kept_dt = dt[, ..remaining],
    drop_log = drop_dt,
    final_vif = final_vif
  )
}
###
# B-4.2.8 run
###
feature_meta <- make_feature_metadata(feature_cols)
metrics_dt <- calc_feature_metrics(x_nzv)
#metrics_dt <- merge(metrics_dt, feature_meta, by = "feature", all.x = TRUE)
head(metrics_dt[order(-normality_score, -mad)])
## run family cap
fam_res <- family_cap_prune(
  dt = x_df,
  metrics_dt = metrics_dt,
  feature_meta = feature_meta,
  group_col = "family_cap_group",
  variability_metric = "mad",
  preference = "normality_then_variability"
)

# 더 넓게 inter_feature_family()기준으로 cap하고 싶다면...
# fam_res_broad <- family_cap_prune(
#   dt = x_imp,
#   metrics_dt = metrics_dt,
#   feature_meta = feature_meta,
#   group_col = "family",
#   variability_metric = "mad",
#   preference = "normality_then_variability"
# )

x_after_family <- fam_res$kept_dt
length(fam_res$kept_features) #360
head(fam_res$drop_log)

## pairwise correlation pruning
metrics_after_family <- metrics_dt[feature %in% colnames(x_after_family)]

cor_res <- correlation_prune(
  dt = x_after_family,
  metrics_dt = metrics_after_family,
  threshold = 0.8,
  method = "pearson",
  variability_metric = "mad",
  preference = "normality_then_variability"
)

x_after_cor <- cor_res$kept_dt
length(cor_res$kept_features) #140
head(cor_res$drop_log)

## VIF pruning
metrics_after_cor <- metrics_dt[feature %in% colnames(x_after_cor)]

vif_res <- vif_prune(
  dt = x_after_cor,
  metrics_dt = metrics_after_cor,
  vif_threshold = 10,
  variability_metric = "mad",
  preference = "normality_then_variability",
  verbose = TRUE
)

x_final <- vif_res$kept_dt
length(vif_res$kept_features) #42
head(vif_res$final_vif[order(-vif)])
head(vif_res$drop_log)

## summary
cat("Original features:     ", ncol(x_df), "\n")
cat("After family cap:      ", ncol(x_after_family), "\n")
cat("After cor pruning:     ", ncol(x_after_cor), "\n")
cat("After VIF pruning:     ", ncol(x_final), "\n")

## save log
dir.create("data/20260309_pilot/results/version2/feature_pruning_logs", showWarnings = FALSE)

fwrite(metrics_dt, "data/20260309_pilot/results/version2/feature_pruning_logs/feature_metrics.csv")

if (nrow(fam_res$drop_log) > 0) {
  fwrite(fam_res$drop_log, "data/20260309_pilot/results/version2/feature_pruning_logs/family_cap_drop_log.csv")
}

if (nrow(cor_res$drop_log) > 0) {
  fwrite(cor_res$drop_log, "data/20260309_pilot/results/version2/feature_pruning_logs/correlation_prune_drop_log.csv")
}

if (nrow(vif_res$drop_log) > 0) {
  fwrite(vif_res$drop_log, "data/20260309_pilot/results/version2/feature_pruning_logs/vif_prune_drop_log.csv")
}

fwrite(data.table(feature = colnames(x_final)),
       "data/20260309_pilot/results/version2/feature_pruning_logs/final_kept_features.csv")



####
# B-5. choose variable feature 
# candidate: MAD, (optional: variance)
####

mad_dt <- data.table(
  feature = colnames(x_final), #x_nzv -> x_final
  mad = sapply(x_final, mad, na.rm = TRUE)
)[order(-mad)]

head(mad_dt, 20)

n_top <- min(20, nrow(mad_dt))
top_features_mad <- mad_dt$feature[1:n_top]

x_top <- x_nzv[, ..top_features_mad]
dim(x_top)

#####
# B-4.3. use clinical data and run FAMD instead of PCA
# B-4.3.1 add clinical data
# skip scaling!
#####

x_top_dt <- data.table(
  row_id = row_id,
  sample = sample_id,
  x_top
)
x_top_dt <- x_top_dt[,-1]
continuous_clinical <- c("Age_IO")
binary_clinical <- c("Sex M: 1  F: 2", "ECOG _PS", "EGFR_Mutation_Status_mod","IO_Line","Previous_palliative_chemo","Previous_palliative_target","PD-L1_TPS_mod")

use_cols <- c('sample',continuous_clinical, binary_clinical)
clinic_df<- merged_dt_anno[, ..use_cols]

mix_df <- merge(clinic_df, x_top_dt,by='sample')

# # molecular continuous
# for (col in molecular_features) {
#   mix_df[[col]] <- as.numeric(x_final[[col]])
# }

# continuous clinical
for (col in continuous_clinical) {
  mix_df[[col]] <- as.numeric(mix_df[[col]])
  mix_df[[col]][is.na(mix_df[[col]])] <- median(mix_df[[col]], na.rm = TRUE)
}

#####
# B-4.3.2 FAMD
#####
mix_df_v2 <- mix_df
rownames(mix_df_v2) <- mix_df_v2$sample
# remove sample column
mix_df_v2$sample <- NULL
# change character to factor
char_cols <- names(mix_df_v2)[sapply(mix_df_v2, is.character)]
char_cols

for (col in char_cols) {
  mix_df_v2[[col]] <- as.factor(mix_df_v2[[col]])
}
# check NA or level
num_cols <- names(mix_df_v2)[sapply(mix_df_v2, is.numeric)]
fac_cols <- names(mix_df_v2)[sapply(mix_df_v2, is.factor)]
check_num_dt <- data.table(
  feature = num_cols,
  n_na = sapply(mix_df_v2[, ..num_cols], function(x) sum(is.na(x))),
  n_nan = sapply(mix_df_v2[, ..num_cols], function(x) sum(is.nan(x))),
  n_inf = sapply(mix_df_v2[, ..num_cols], function(x) sum(is.infinite(x))),
  sd = sapply(mix_df_v2[, ..num_cols], function(x) sd(x, na.rm = TRUE))
)
check_fac_dt <- data.table(
  feature = fac_cols,
  n_na = sapply(mix_df_v2[, ..fac_cols], function(x) sum(is.na(x))),
  n_levels = sapply(mix_df_v2[, ..fac_cols], function(x) nlevels(droplevels(x)))
)

mix_df_v2$`PD-L1_TPS_mod` <- as.character(mix_df_v2$`PD-L1_TPS_mod`)
mix_df_v2$`PD-L1_TPS_mod`[is.na(mix_df_v2$`PD-L1_TPS_mod`)] <- "Unknown"
mix_df_v2$`PD-L1_TPS_mod` <- factor(mix_df_v2$`PD-L1_TPS_mod`)

# run
famd_res <- FAMD(mix_df_v2, graph = FALSE)

# elbow와 비슷하게 보려면
eig_dt <- data.table(
  Dim = seq_len(nrow(famd_res$eig)),
  eigenvalue = famd_res$eig[,1],
  pct_var = famd_res$eig[, 2],
  cum_var = famd_res$eig[, 3]
)

eig_dt

library(FactoMineR)
desc <- dimdesc(famd_res, axes = 1:5)

desc$Dim.1$category

ggplot(eig_dt, aes(Dim, pct_var)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  labs(
    x = "FAMD dimension",
    y = "Percent variance explained",
    title = "FAMD elbow plot"
  )#20260428_version2_trial2_FAMD

famd_coord <- famd_res$ind$coord
rownames(famd_coord) <- mix_df$sample

famd_dt <- data.table(
  sample = rownames(famd_coord) ,
  famd_res$ind$coord[, 1:n_dim, drop = FALSE]
)

setnames(
  famd_dt,
  old = names(famd_dt)[-1],
  new = paste0("Dim", 1:n_dim)
)
famd_plot_dt <- merge(
  famd_dt,
  meta,
  by.x = "sample",
  by.y = "sample_id"
)
famd_plot_dt_batch <- merge(famd_plot_dt,old[,c('patient_id','batch')],by.x='sample',by.y='patient_id',all.x=T)
ggplot(famd_plot_dt_batch, aes(x = Dim1, y = Dim2, color = `Binarized response`,label=batch)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3)+
  theme_bw() +
  labs(
    x = "FAMD Dim1",
    y = "FAMD Dim2",
    title = "FAMD plot"
  )#20260428_version2_trial2_FAMD_FD12 #20260428_version2_trial2_FAMD_FD23

#####
# B-4.3.3 component 선택 후 clustering
#####
n_dim <- 5
emb <- famd_res$ind$coord[, 2:n_dim, drop = FALSE]
# neighbor
k_nn <- 10 #5
knn_res <- get.knn(emb, k = k_nn)

# go to ->>>> B-9. k-nearest neighbors graph 

#####
# B-6. scaling
#####
x_scaled <- scale(x_top) 
rownames(x_scaled) <- sample_id
dim(x_scaled)


#####
# B-7. PCA
####
pca_res <- prcomp(x_scaled, center = FALSE, scale. = FALSE)
summary(pca_res)

pca_dt <- data.table(
  sample = rownames(x_scaled),
  pca_res$x[, 1:20, drop = FALSE]
)

setnames(pca_dt, old = paste0("PC", 1:20), new = paste0("PC", 1:20))

ggplot(pca_dt, aes(PC1, PC2, label = sample)) +
  geom_point(size = 3) +
  theme_bw()
pca_dt_plot <- merge(pca_dt, meta, by.x='sample',by.y='sample_id')
pca_dt_batch <- merge(pca_dt_plot, old[,c('patient_id','batch')],by.x='sample',by.y='patient_id')
ggplot(pca_dt_batch, aes(PC1, PC2, label = batch,col=batch)) +
  geom_text(vjust = -0.5, size = 3)+
  geom_point(size = 3) +
  theme_bw()

####
dt <- copy(pca_dt_batch)

pc_cols <- paste0("PC", 1:20)

# ID나 이미 PCA에서 나온 column은 제외
exclude_cols <- c(
  "sample", "Study ID", "Sample ID",
  pc_cols
)

meta_cols <- setdiff(names(dt), exclude_cols)

# 값이 2개 이상 있는 column만 사용
meta_cols <- meta_cols[
  sapply(dt[, ..meta_cols], function(x) data.table::uniqueN(na.omit(x)) >= 2)
]

run_adonis_one <- function(dt, var, pcs = pc_cols, nperm = 999) {
  sub <- dt[, c(pcs, var), with = FALSE]
  sub <- na.omit(sub)
  
  if (nrow(sub) < 5) return(NULL)
  if (data.table::uniqueN(sub[[var]]) < 2) return(NULL)
  
  # PC는 numeric 보장
  for (p in pcs) sub[[p]] <- as.numeric(sub[[p]])
  
  # metadata column을 안전한 이름으로 복사
  if (is.character(sub[[var]]) || is.factor(sub[[var]])) {
    sub[, meta_var := factor(get(var))]
  } else {
    sub[, meta_var := as.numeric(get(var))]
  }
  
  X <- as.matrix(sub[, ..pcs])
  d <- dist(X)
  
  fit <- vegan::adonis2(
    d ~ meta_var,
    data = as.data.frame(sub),
    permutations = nperm
  )
  
  data.table(
    var = var,
    pc_set = paste(pcs, collapse = "+"),
    n = nrow(sub),
    n_level = data.table::uniqueN(sub$meta_var),
    R2 = fit$R2[1],
    F = fit$F[1],
    p = fit$`Pr(>F)`[1]
  )
}

global_pc_res <- rbindlist(
  lapply(meta_cols, function(v) run_adonis_one(dt, v, pcs = pc_cols)),
  fill = TRUE
)

global_pc_res[, p_adj := p.adjust(p, method = "BH")]
setorder(global_pc_res, -R2)

global_pc_res
global_pc_res[order(-R2)]
## 2.어떤 PC 하나가 어떤 metadata를 잘 설명하는지 heatmap
single_pc_res <- rbindlist(
  lapply(meta_cols, function(v) {
    rbindlist(
      lapply(pc_cols, function(pc) {
        run_adonis_one(dt, v, pcs = pc)
      }),
      fill = TRUE
    )
  }),
  fill = TRUE
)

single_pc_res[, PC := pc_set]
single_pc_res[, p_adj := p.adjust(p, method = "BH")]
single_pc_res[, neglog10_padj := -log10(p_adj)]

single_pc_res[order(-R2)]
x_order <- paste0('PC',seq(1,20))
ggplot(single_pc_res, aes(x = factor(PC,levels=x_order), y = reorder(var, R2), fill = R2)) +
  geom_tile() +
  theme_bw() +
  labs(
    x = "PC",
    y = "Metadata column",
    fill = "R2",
    title = "How strongly each PC is associated with each metadata column"
  )+theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))#5x6 #20260428_version2_trial2 #20260428_version2_trial2_exceptBatch
####

# find features explaining specific PC
library(pheatmap)
library(data.table)
library(pheatmap)

library(data.table)
library(pheatmap)

make_annotation_colors <- function(annotation_row = NULL, annotation_col = NULL) {
  ann_colors <- list()
  
  # row annotation: direction 색 명시
  if (!is.null(annotation_row) && "direction" %in% names(annotation_row)) {
    ann_colors$direction <- c(
      positive = "#D73027",
      negative = "#4575B4"
    )
  }
  
  # row annotation: loading 색 명시
  # positive loading = red, negative loading = blue
  if (!is.null(annotation_row) && "loading" %in% names(annotation_row)) {
    ann_colors$loading <- grDevices::colorRampPalette(
      c("#4575B4", "white", "#D73027")
    )(100)
  }
  
  add_discrete_colors <- function(df) {
    if (is.null(df)) return(NULL)
    
    for (nm in names(df)) {
      x <- df[[nm]]
      
      # 이미 지정한 annotation은 skip
      if (nm %in% names(ann_colors)) next
      
      # numeric은 continuous로 처리하게 둠
      if (is.numeric(x)) next
      
      vals <- sort(unique(as.character(na.omit(x))))
      if (length(vals) == 0) next
      
      cols <- grDevices::hcl.colors(
        n = length(vals),
        palette = "Dynamic"
      )
      names(cols) <- vals
      
      ann_colors[[nm]] <<- cols
    }
  }
  
  add_discrete_colors(annotation_row)
  add_discrete_colors(annotation_col)
  
  ann_colors
}


plot_pc_feature_heatmap <- function(
    pca_res,
    x_input,
    pc = "PC4",
    top_n = 5,
    meta_dt = NULL,
    sample_id_col = "sample",
    annotation_cols = NULL,
    scale_by_feature = FALSE
) {
  if (is.null(rownames(x_input))) {
    stop("x_input must have rownames as sample IDs.")
  }
  if (is.null(colnames(x_input))) {
    stop("x_input must have colnames as feature names.")
  }
  if (!(pc %in% colnames(pca_res$rotation))) {
    stop(sprintf("'%s' not found in pca_res$rotation.", pc))
  }
  
  # loading table
  loadings_dt <- as.data.table(
    pca_res$rotation,
    keep.rownames = "feature"
  )
  
  pc_dt <- loadings_dt[, .(
    feature,
    loading = get(pc)
  )]
  
  # top positive / negative feature
  top_pos <- pc_dt[order(-loading)][1:min(top_n, .N)]
  top_pos[, direction := "positive"]
  top_pos[, rank := seq_len(.N)]
  
  top_neg <- pc_dt[order(loading)][1:min(top_n, .N)]
  top_neg[, direction := "negative"]
  top_neg[, rank := seq_len(.N)]
  
  selected_dt <- rbind(top_pos, top_neg, fill = TRUE)
  
  feature_order <- c(top_pos$feature, top_neg$feature)
  feature_order <- unique(feature_order)
  
  # sample order by PC score
  sample_scores <- data.table(
    sample = rownames(pca_res$x),
    PC_score = pca_res$x[, pc]
  )
  setorder(sample_scores, -PC_score)
  sample_order <- sample_scores$sample
  
  common_samples <- intersect(sample_order, rownames(x_input))
  common_features <- intersect(feature_order, colnames(x_input))
  
  if (length(common_samples) == 0) {
    stop("No overlapping samples between pca_res$x and x_input.")
  }
  if (length(common_features) == 0) {
    stop("No selected features found in x_input columns.")
  }
  
  mat <- x_input[common_samples, common_features, drop = FALSE]
  mat <- mat[, feature_order[feature_order %in% colnames(mat)], drop = FALSE]
  
  # rows = features, columns = samples
  heatmap_mat <- t(mat)
  
  if (scale_by_feature) {
    heatmap_mat <- t(scale(t(heatmap_mat)))
    heatmap_mat[is.na(heatmap_mat)] <- 0
  }
  
  # row annotation
  row_annot <- selected_dt[match(rownames(heatmap_mat), feature), .(
    direction,
    loading
  )]
  row_annot <- as.data.frame(row_annot)
  rownames(row_annot) <- rownames(heatmap_mat)
  row_annot$direction <- factor(row_annot$direction, levels = c("positive", "negative"))
  
  # column annotation
  if (!is.null(meta_dt)) {
    meta2 <- as.data.table(meta_dt)
    
    if (!(sample_id_col %in% names(meta2))) {
      stop(sprintf("sample_id_col '%s' not found in meta_dt.", sample_id_col))
    }
    
    meta2 <- unique(meta2, by = sample_id_col)
    meta2 <- meta2[get(sample_id_col) %in% colnames(heatmap_mat)]
    
    if (!is.null(annotation_cols)) {
      keep_cols <- c(sample_id_col, annotation_cols[annotation_cols %in% names(meta2)])
    } else {
      keep_cols <- sample_id_col
    }
    
    meta2 <- meta2[, ..keep_cols]
    
    meta2 <- merge(
      meta2,
      sample_scores,
      by.x = sample_id_col,
      by.y = "sample",
      all.x = TRUE
    )
    
    col_annot <- as.data.frame(meta2)
    rownames(col_annot) <- col_annot[[sample_id_col]]
    col_annot[[sample_id_col]] <- NULL
    
    col_annot <- col_annot[colnames(heatmap_mat), , drop = FALSE]
    
    # character columns를 factor로 명확히 변환
    for (nm in names(col_annot)) {
      if (is.character(col_annot[[nm]])) {
        col_annot[[nm]] <- factor(col_annot[[nm]])
      }
    }
  } else {
    col_annot <- data.frame(
      PC_score = sample_scores[match(colnames(heatmap_mat), sample), PC_score]
    )
    rownames(col_annot) <- colnames(heatmap_mat)
  }
  
  ann_colors <- make_annotation_colors(
    annotation_row = row_annot,
    annotation_col = col_annot
  )
  
  pheatmap(
    mat = heatmap_mat,
    annotation_col = col_annot,
    annotation_row = row_annot,
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    annotation_legend = TRUE,
    legend = TRUE,
    main = paste0(
      pc, ": top +", top_n, " / -", top_n, " feature heatmap"
    )
  )
  invisible(list(
    selected_features = selected_dt,
    heatmap_mat = heatmap_mat,
    sample_scores = sample_scores,
    annotation_col = col_annot,
    annotation_row = row_annot
  ))
}
res_pc4 <- plot_pc_feature_heatmap(
  pca_res = pca_res,
  x_input = x_scaled,
  pc = "PC4",
  top_n = 5,
  meta_dt = pca_plot_dt_batch,
  sample_id_col = "sample",
  annotation_cols = c("cluster", "Binarized response","Histology_mod"),
  scale_by_feature = FALSE
)#3x12 #20260428_version2_trial2_PC4 #20260428_version2_trial2_PC3

####

var_explained <- pca_res$sdev^2
pct_var <- var_explained / sum(var_explained) * 100
cum_var <- cumsum(pct_var)

pca_var_dt <- data.table(
  PC = seq_along(pct_var),
  pct_var = pct_var,
  cum_var = cum_var
)

ggplot(pca_var_dt, aes(x = PC, y = pct_var)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  labs(
    x = "Principal component (PC)",
    y = "Percent variance explained",
    title = "PCA elbow plot"
  )#3x4 #20260428_version2_trial2_pc_variance_explained

npc <-seq(2,11)
npc <- npc[npc!=3]
pc_use <- pca_res$x[, npc, drop = FALSE]
dim(pc_use)


#####
# B-8. neighbors
#####
#k_nn <- min(10, nrow(pc_use) - 1)
k_nn <- 10
knn_res <- get.knn(pc_use, k = k_nn)
str(knn_res)

#####
# B-9. k-nearest neighbors graph 
#####
edge_df <- data.table()

for (i in seq_len(nrow(knn_res$nn.index))) {
  tmp <- data.table(
    from = i,
    to = knn_res$nn.index[i, ],
    weight = 1 / (1 + knn_res$nn.dist[i, ])
  )
  edge_df <- rbind(edge_df, tmp)
}

head(edge_df)

# scaling -> pca
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = seq_len(nrow(pc_use))))
# famb
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = seq_len(nrow(emb))))
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "mean")

g
#####
# B-10. Louvain clustering
#####
clu <- cluster_louvain(g, weights = E(g)$weight, resolution =0.5) # resolution
membership_vec <- membership(clu)

table(membership_vec)

# scaling -> pca
cluster_dt <- data.table(
  sample = rownames(x_scaled),
  cluster = factor(membership_vec)
)
# famb
cluster_dt <- data.table(
  sample = mix_df$sample,
  cluster = factor(membership_vec)
)


head(cluster_dt)

# scaling -> pca
pca_plot_dt <- merge(pca_dt_batch, cluster_dt, by = "sample", all.x = TRUE)

ggplot(pca_plot_dt, aes(PC1, PC2, color = cluster, shape=`Binarized response`,label = sample)) +
  geom_point(size = 3) +
  theme_bw()


#####
# B-11. umap
#####
set.seed(123)
# scaling -> pca
umap_res <- umap(
  pc_use,
  n_neighbors = min(10, nrow(pc_use) - 1),
  min_dist = 0.3,
  metric = "euclidean"
)
umap_dt <- data.table(
  sample = rownames(x_scaled),
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2]
)
# famb
umap_res <- umap(
  emb,
  n_neighbors = min(10, nrow(emb) - 1),
  min_dist = 0.3,
  metric = "euclidean"
)
umap_dt <- data.table(
  sample = mix_df$sample,
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2]
)

umap_plot_dt <- merge(umap_dt, cluster_dt, by = "sample", all.x = TRUE)
head(umap_plot_dt)

table(umap_plot_dt_meta$cluster,umap_plot_dt_meta$`Binarized response`)
old <- fread('data/20260309_pilot/results/version1/subtype_discovery/base_run_v4/patient_features.csv')
old$patient_id <- gsub('_','-',old$patient_id)
umap_plot_dt_batch <- merge(umap_plot_dt_meta, old[,c('patient_id','batch')],by.x='sample',by.y='patient_id')
ggplot(umap_plot_dt_batch, aes(UMAP1, UMAP2, color = cluster,label = batch)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw()#20260428_version2_trial2_umap_v2 #20260428_version2_trial2_FAMD_umap #20260428_version2_trial2_FAMD_umap_v2




#####
# A-FAMD version
####
#####
# B-4.3.2 FAMD
#####
x_mat_v2 <- x_mat
x_mat_v2 <- as.data.frame(x_mat_v2)
x_mat_v2$sample <- rownames(x_mat_v2)

continuous_clinical <- c("Age_IO")
binary_clinical <- c("Sex M: 1  F: 2", "ECOG _PS", "EGFR_Mutation_Status_mod","IO_Line","Previous_palliative_chemo","Previous_palliative_target","PD-L1_TPS_mod")

use_cols <- c('sample',continuous_clinical, binary_clinical)
clinic_df<- merged_dt_anno[, ..use_cols]

mix_df_A_famd <- merge(clinic_df, x_mat_v2,by='sample')

# # molecular continuous
# for (col in molecular_features) {
#   mix_df[[col]] <- as.numeric(x_final[[col]])
# }

# continuous clinical
for (col in continuous_clinical) {
  mix_df_A_famd[[col]] <- as.numeric(mix_df_A_famd[[col]])
  mix_df_A_famd[[col]][is.na(mix_df_A_famd[[col]])] <- median(mix_df_A_famd[[col]], na.rm = TRUE)
}


mix_df_A_famd_v2 <- mix_df_A_famd
rownames(mix_df_A_famd_v2) <- mix_df_A_famd_v2$sample
# remove sample column
mix_df_A_famd_v2$sample <- NULL
# change character to factor
char_cols <- names(mix_df_A_famd_v2)[sapply(mix_df_A_famd_v2, is.character)]
char_cols

for (col in char_cols) {
  mix_df_A_famd_v2[[col]] <- as.factor(mix_df_A_famd_v2[[col]])
}
# check NA or level
num_cols <- names(mix_df_A_famd_v2)[sapply(mix_df_A_famd_v2, is.numeric)]
fac_cols <- names(mix_df_A_famd_v2)[sapply(mix_df_A_famd_v2, is.factor)]
check_num_dt <- data.table(
  feature = num_cols,
  n_na = sapply(mix_df_A_famd_v2[, ..num_cols], function(x) sum(is.na(x))),
  n_nan = sapply(mix_df_A_famd_v2[, ..num_cols], function(x) sum(is.nan(x))),
  n_inf = sapply(mix_df_A_famd_v2[, ..num_cols], function(x) sum(is.infinite(x))),
  sd = sapply(mix_df_A_famd_v2[, ..num_cols], function(x) sd(x, na.rm = TRUE))
)
check_fac_dt <- data.table(
  feature = fac_cols,
  n_na = sapply(mix_df_A_famd_v2[, ..fac_cols], function(x) sum(is.na(x))),
  n_levels = sapply(mix_df_A_famd_v2[, ..fac_cols], function(x) nlevels(droplevels(x)))
)

mix_df_A_famd_v2$`PD-L1_TPS_mod` <- as.character(mix_df_A_famd_v2$`PD-L1_TPS_mod`)
mix_df_A_famd_v2$`PD-L1_TPS_mod`[is.na(mix_df_A_famd_v2$`PD-L1_TPS_mod`)] <- "Unknown"
mix_df_A_famd_v2$`PD-L1_TPS_mod` <- factor(mix_df_A_famd_v2$`PD-L1_TPS_mod`)

# run
x_famd <- as.data.frame(mix_df_A_famd_v2)
rownames(x_famd) <- mix_df_A_famd$sample
famd_res <- FAMD(x_famd, graph = FALSE)

# elbow와 비슷하게 보려면
eig_dt <- data.table(
  Dim = seq_len(nrow(famd_res$eig)),
  eigenvalue = famd_res$eig[,1],
  pct_var = famd_res$eig[, 2],
  cum_var = famd_res$eig[, 3]
)

eig_dt

library(FactoMineR)
desc <- dimdesc(famd_res, axes = 1:5)


ggplot(eig_dt, aes(Dim, pct_var)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  labs(
    x = "FAMD dimension",
    y = "Percent variance explained",
    title = "FAMD elbow plot"
  )#3x4 #20260429_version2_trial1_FAMD_elbow

famd_coord <- famd_res$ind$coord
rownames(famd_coord) <- mix_df$sample

famd_dt <- data.table(
  sample = rownames(famd_coord) ,
  famd_res$ind$coord[, 1:n_dim, drop = FALSE]
)

setnames(
  famd_dt,
  old = names(famd_dt)[-1],
  new = paste0("Dim", 1:n_dim)
)
famd_plot_dt <- merge(
  famd_dt,
  meta,
  by.x = "sample",
  by.y = "sample_id"
)
famd_plot_dt_batch <- merge(famd_plot_dt,old[,c('patient_id','batch')],by.x='sample',by.y='patient_id',all.x=T)
ggplot(famd_plot_dt_batch, aes(x = Dim2, y = Dim3, color = `Binarized response`,label=batch)) +
  geom_point(size = 3,alpha=0.6) +
  geom_text(vjust = -0.5, size = 2)+
  theme_bw() +
  labs(
    x = "FAMD Dim2",
    y = "FAMD Dim3",
    title = "FAMD plot"
  )#20260429_version2_trial1_FAMD_FD12 #20260429_version2_trial1_FAMD_FD23

# visualziatino
make_annotation_colors <- function(annotation_row = NULL, annotation_col = NULL) {
  ann_colors <- list()
  
  if (!is.null(annotation_row)) {
    if ("direction" %in% names(annotation_row)) {
      ann_colors$direction <- c(
        positive = "#D73027",
        negative = "#4575B4"
      )
    }
    if ("item_type" %in% names(annotation_row)) {
      vals <- sort(unique(as.character(na.omit(annotation_row$item_type))))
      cols <- grDevices::hcl.colors(length(vals), palette = "Set 2")
      names(cols) <- vals
      ann_colors$item_type <- cols
    }
  }
  
  add_discrete_colors <- function(df) {
    if (is.null(df)) return(NULL)
    
    for (nm in names(df)) {
      x <- df[[nm]]
      if (nm %in% names(ann_colors)) next
      if (is.numeric(x)) next
      
      vals <- sort(unique(as.character(na.omit(x))))
      if (length(vals) == 0) next
      
      cols <- grDevices::hcl.colors(
        n = length(vals),
        palette = "Dynamic"
      )
      names(cols) <- vals
      
      ann_colors[[nm]] <<- cols
    }
  }
  
  add_discrete_colors(annotation_row)
  add_discrete_colors(annotation_col)
  
  ann_colors
}
get_score_col <- function(dt) {
  cand <- c("Estimate", "estimate", "correlation", "coord", "v.test")
  hit <- cand[cand %in% names(dt)]
  if (length(hit) > 0) return(hit[1])
  
  num_cols <- names(dt)[sapply(dt, is.numeric)]
  num_cols <- setdiff(num_cols, "p.value")
  if (length(num_cols) == 0) {
    stop("Could not find a numeric score column in dimdesc table.")
  }
  num_cols[1]
}
parse_modality_label <- function(label, input_colnames = NULL) {
  label <- as.character(label)
  
  # 가장 흔한 형식
  if (grepl("=", label, fixed = TRUE)) {
    sp <- strsplit(label, "=", fixed = TRUE)[[1]]
    return(list(feature = sp[1], modality = paste(sp[-1], collapse = "=")))
  }
  if (grepl(":", label, fixed = TRUE)) {
    sp <- strsplit(label, ":", fixed = TRUE)[[1]]
    return(list(feature = sp[1], modality = paste(sp[-1], collapse = ":")))
  }
  
  # colname prefix 기반 추정
  if (!is.null(input_colnames)) {
    hits <- input_colnames[startsWith(label, paste0(input_colnames, "_"))]
    if (length(hits) > 0) {
      feat <- hits[which.max(nchar(hits))]
      mod <- sub(paste0("^", feat, "_"), "", label)
      return(list(feature = feat, modality = mod))
    }
  }
  
  # fallback
  list(feature = label, modality = NA_character_)
}
extract_famd_dim_features <- function(
    desc_axis,
    x_input,
    top_n = 5,
    include_quanti = TRUE,
    include_category = TRUE
) {
  out_list <- list()
  
  # quantitative variables
  if (include_quanti && "quanti" %in% names(desc_axis) && !is.null(desc_axis$quanti)) {
    qdt <- as.data.table(desc_axis$quanti, keep.rownames = "item")
    if (nrow(qdt) > 0) {
      score_col <- get_score_col(qdt)
      qdt[, score := get(score_col)]
      qdt[, pval := if ("p.value" %in% names(qdt)) `p.value` else NA_real_]
      qdt[, item_type := "quanti"]
      qdt[, feature := item]
      qdt[, modality := NA_character_]
      out_list[["quanti"]] <- qdt[, .(item, feature, modality, item_type, score, pval)]
    }
  }
  
  # category / modality-level descriptors
  cat_name <- intersect(c("category", "categories"), names(desc_axis))
  if (include_category && length(cat_name) > 0) {
    cdt <- as.data.table(desc_axis[[cat_name[1]]], keep.rownames = "item")
    if (nrow(cdt) > 0) {
      score_col <- get_score_col(cdt)
      cdt[, score := get(score_col)]
      cdt[, pval := if ("p.value" %in% names(cdt)) `p.value` else NA_real_]
      cdt[, item_type := "category"]
      
      parsed <- lapply(cdt$item, parse_modality_label, input_colnames = colnames(x_input))
      cdt[, feature := vapply(parsed, `[[`, character(1), "feature")]
      cdt[, modality := vapply(parsed, `[[`, character(1), "modality")]
      
      out_list[["category"]] <- cdt[, .(item, feature, modality, item_type, score, pval)]
    }
  }
  
  if (length(out_list) == 0) {
    stop("No usable descriptors found in dimdesc for this axis.")
  }
  
  all_dt <- rbindlist(out_list, fill = TRUE)
  
  # top positive / negative
  top_pos <- all_dt[order(-score)][1:min(top_n, .N)]
  top_pos[, direction := "positive"]
  top_pos[, rank_within_direction := seq_len(.N)]
  
  top_neg <- all_dt[order(score)][1:min(top_n, .N)]
  top_neg[, direction := "negative"]
  top_neg[, rank_within_direction := seq_len(.N)]
  
  selected_dt <- rbind(top_pos, top_neg, fill = TRUE)
  
  # 중복 item 제거
  selected_dt <- unique(selected_dt, by = "item")
  
  # positive 먼저, negative 나중
  pos_items <- top_pos$item
  neg_items <- top_neg$item
  item_order <- unique(c(pos_items, neg_items))
  selected_dt <- selected_dt[match(item_order, item)]
  
  selected_dt
}

plot_famd_feature_heatmap <- function(
    famd_res,
    x_input,
    desc = NULL,
    dim = 1,
    top_n = 5,
    meta_dt = NULL,
    sample_id_col = "sample",
    annotation_cols = NULL,
    include_quanti = TRUE,
    include_category = TRUE,
    scale_quanti_rows = TRUE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE
) {
  if (is.null(rownames(x_input))) {
    stop("x_input must have rownames as sample IDs.")
  }
  
  if (is.null(desc)) {
    desc <- dimdesc(famd_res, axes = dim)
  }
  
  dim_name <- paste0("Dim.", dim)
  
  # dimdesc axis object 가져오기
  desc_axis <- NULL
  if (dim_name %in% names(desc)) {
    desc_axis <- desc[[dim_name]]
  } else if (length(desc) >= dim) {
    desc_axis <- desc[[dim]]
  } else {
    stop(sprintf("Could not find %s in dimdesc result.", dim_name))
  }
  
  # sample scores
  if (is.null(famd_res$ind$coord) || !(dim_name %in% colnames(famd_res$ind$coord))) {
    stop(sprintf("'%s' not found in famd_res$ind$coord.", dim_name))
  }
  
  sample_scores <- data.table(
    sample = rownames(famd_res$ind$coord),
    Dim_score = famd_res$ind$coord[, dim_name]
  )
  setorder(sample_scores, -Dim_score)
  sample_order <- sample_scores$sample
  
  # selected features/categories
  selected_dt <- extract_famd_dim_features(
    desc_axis = desc_axis,
    x_input = x_input,
    top_n = top_n,
    include_quanti = include_quanti,
    include_category = include_category
  )
  
  common_samples <- intersect(sample_order, rownames(x_input))
  if (length(common_samples) == 0) {
    stop("No overlapping samples between FAMD result and x_input.")
  }
  
  # feature x sample matrix 만들기
  row_list <- list()
  
  for (i in seq_len(nrow(selected_dt))) {
    feature_i <- selected_dt$feature[i]
    modality_i <- selected_dt$modality[i]
    type_i <- selected_dt$item_type[i]
    item_i <- selected_dt$item[i]
    
    if (!(feature_i %in% colnames(x_input))) next
    
    vec <- x_input[common_samples, feature_i]
    
    if (type_i == "quanti") {
      val <- as.numeric(vec)
    } else if (type_i == "category") {
      val <- as.numeric(as.character(vec) == modality_i)
    } else {
      next
    }
    
    row_list[[item_i]] <- val
  }
  
  if (length(row_list) == 0) {
    stop("No selected rows could be built for the heatmap.")
  }
  
  heatmap_mat <- do.call(rbind, row_list)
  colnames(heatmap_mat) <- common_samples
  
  # quantitative row만 scale
  row_annot <- selected_dt[match(rownames(heatmap_mat), item), .(
    feature,
    modality,
    item_type,
    direction,
    score,
    pval
  )]
  row_annot <- as.data.frame(row_annot)
  rownames(row_annot) <- rownames(heatmap_mat)
  row_annot$direction <- factor(row_annot$direction, levels = c("positive", "negative"))
  row_annot$item_type <- factor(row_annot$item_type)
  
  if (scale_quanti_rows) {
    quanti_rows <- which(row_annot$item_type == "quanti")
    if (length(quanti_rows) > 0) {
      heatmap_mat[quanti_rows, ] <- t(scale(t(heatmap_mat[quanti_rows, , drop = FALSE])))
      heatmap_mat[is.na(heatmap_mat)] <- 0
    }
  }
  
  # column annotation
  if (!is.null(meta_dt)) {
    meta2 <- as.data.table(meta_dt)
    
    if (!(sample_id_col %in% names(meta2))) {
      stop(sprintf("sample_id_col '%s' not found in meta_dt.", sample_id_col))
    }
    
    meta2 <- unique(meta2, by = sample_id_col)
    meta2 <- meta2[get(sample_id_col) %in% colnames(heatmap_mat)]
    
    if (!is.null(annotation_cols)) {
      keep_cols <- c(sample_id_col, annotation_cols[annotation_cols %in% names(meta2)])
    } else {
      keep_cols <- sample_id_col
    }
    
    meta2 <- meta2[, ..keep_cols]
    
    meta2 <- merge(
      meta2,
      sample_scores,
      by.x = sample_id_col,
      by.y = "sample",
      all.x = TRUE
    )
    
    col_annot <- as.data.frame(meta2)
    rownames(col_annot) <- col_annot[[sample_id_col]]
    col_annot[[sample_id_col]] <- NULL
    col_annot <- col_annot[colnames(heatmap_mat), , drop = FALSE]
    
    for (nm in names(col_annot)) {
      if (is.character(col_annot[[nm]])) {
        col_annot[[nm]] <- factor(col_annot[[nm]])
      }
    }
  } else {
    col_annot <- data.frame(
      Dim_score = sample_scores[match(colnames(heatmap_mat), sample), Dim_score]
    )
    rownames(col_annot) <- colnames(heatmap_mat)
  }
  
  ann_colors <- make_annotation_colors(
    annotation_row = row_annot[, c("item_type", "direction"), drop = FALSE],
    annotation_col = col_annot
  )
  
  pheatmap(
    mat = heatmap_mat,
    annotation_col = col_annot,
    annotation_row = row_annot[, c("item_type", "direction"), drop = FALSE],
    annotation_colors = ann_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_colnames = show_colnames,
    annotation_legend = TRUE,
    legend = TRUE,
    border_color = NA,
    main = paste0(
      dim_name,
      ": top +", top_n,
      " / -", top_n,
      " FAMD descriptors"
    )
  )
  
  invisible(list(
    selected_features = selected_dt,
    heatmap_mat = heatmap_mat,
    sample_scores = sample_scores,
    annotation_col = col_annot,
    annotation_row = row_annot
  ))
}

plot_famd_feature_heatmaps <- function(
    famd_res,
    x_input,
    dims = 1:5,
    desc = NULL,
    top_n = 5,
    meta_dt = NULL,
    sample_id_col = "sample",
    annotation_cols = NULL,
    include_quanti = TRUE,
    include_category = TRUE,
    scale_quanti_rows = TRUE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE
) {
  if (is.null(desc)) {
    desc <- dimdesc(famd_res, axes = dims)
  }
  
  res_list <- vector("list", length(dims))
  names(res_list) <- paste0("Dim.", dims)
  
  for (i in seq_along(dims)) {
    d <- dims[i]
    res_list[[i]] <- plot_famd_feature_heatmap(
      famd_res = famd_res,
      x_input = x_input,
      desc = desc,
      dim = d,
      top_n = top_n,
      meta_dt = meta_dt,
      sample_id_col = sample_id_col,
      annotation_cols = annotation_cols,
      include_quanti = include_quanti,
      include_category = include_category,
      scale_quanti_rows = scale_quanti_rows,
      show_colnames = show_colnames,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols
    )
  }
  
  invisible(res_list)
}
meta$PFS_6mo <- pmin(meta$`PFS (Days)`,180)
res_dim1 <- plot_famd_feature_heatmap(
  famd_res = famd_res,
  x_input = x_famd,
  desc = desc,
  dim = 1,
  top_n = 5,
  meta_dt=meta,
  sample_id_col = 'sample_id',
  annotation_cols = c('Binarized response','RECIST response','PD-L1_TPS_mod','PFS_6mo')
)#4x12 #20260429_version2_trial1_FAMD_FD2_explain #20260429_version2_trial1_FAMD_FD1_explain
#####
# B-4.3.3 component 선택 후 clustering
#####
n_dim <- 5
emb <- famd_res$ind$coord[, 2:n_dim, drop = FALSE]
# neighbor
k_nn <- 5
knn_res <- get.knn(emb, k = k_nn)
#####
# B-9. k-nearest neighbors graph 
#####
edge_df <- data.table()

for (i in seq_len(nrow(knn_res$nn.index))) {
  tmp <- data.table(
    from = i,
    to = knn_res$nn.index[i, ],
    weight = 1 / (1 + knn_res$nn.dist[i, ])
  )
  edge_df <- rbind(edge_df, tmp)
}

head(edge_df)

# famb
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = seq_len(nrow(emb))))
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "mean")

g
#####
# B-10. Louvain clustering
#####
clu <- cluster_louvain(g, weights = E(g)$weight, resolution =1) # resolution
membership_vec <- membership(clu)

table(membership_vec)

# famb
cluster_dt <- data.table(
  sample = mix_df_A_famd$sample,
  cluster = factor(membership_vec)
)


head(cluster_dt)


#####
# B-11. umap
#####
set.seed(123)
# famb
umap_res <- umap(
  emb,
  n_neighbors = min(10, nrow(emb) - 1),
  min_dist = 0.3,
  metric = "euclidean"
)
umap_dt <- data.table(
  sample = mix_df_A_famd$sample,
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2]
)

umap_plot_dt <- merge(umap_dt, cluster_dt, by = "sample", all.x = TRUE)
head(umap_plot_dt)
umap_plot_dt_meta <-merge(umap_plot_dt, meta, by.x='sample',by.y='sample_id',all.x=T)

umap_plot_dt_batch <- merge(umap_plot_dt_meta, old[,c('patient_id','batch')],by.x='sample',by.y='patient_id')
ggplot(umap_plot_dt_batch, aes(UMAP1, UMAP2, color = cluster,label = batch)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw()#20260429_version2_trial1_FAMD_umap
ggplot(umap_plot_dt_batch, aes(UMAP1, UMAP2, color = `Binarized response`,label = batch,shape=cluster)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw()#20260429_version2_trial1_FAMD_umap_binarized_response
umap_plot_dt_batch$PFS_180 <- pmin(umap_plot_dt_batch$`PFS (Days)`, 180)
ggplot(umap_plot_dt_batch, aes(UMAP1, UMAP2, color = PFS_180,label = batch,shape=cluster)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw()#20260429_version2_trial1_FAMD_umap_PFS_6mo
ggplot(umap_plot_dt_batch, aes(UMAP1, UMAP2, color = `RECIST response`,label = batch,shape=cluster)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw()#20260429_version2_trial1_FAMD_umap_RECIST
ggplot(umap_plot_dt_batch, aes(UMAP1, UMAP2, color = `PD-L1_TPS_mod`,label = batch,shape=cluster)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +scale_color_manual(values=c('<1'='yellow','1-49'='orange','≥50'='red','NA'='lightgrey'))+
  theme_bw()#20260429_version2_trial1_FAMD_umap_PDL1_TPS


# cluster 별 marker
mix_dt <- as.data.table(mix_df_A_famd)

# sample column이 있어야 함
stopifnot("sample" %in% names(mix_dt))
stopifnot("sample" %in% names(cluster_dt))
stopifnot("cluster" %in% names(cluster_dt))

analysis_dt <- merge(
  cluster_dt,
  mix_dt,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

analysis_dt[, cluster := factor(cluster)]

table(analysis_dt$cluster)

head(membership_vec)
head(mix_df_A_famd$sample)
length(membership_vec)
nrow(mix_df_A_famd)

id_cols <- c("sample", "cluster")

feature_cols <- setdiff(names(analysis_dt), id_cols)

numeric_cols <- feature_cols[sapply(analysis_dt[, ..feature_cols], is.numeric)]

categorical_cols <- feature_cols[
  sapply(analysis_dt[, ..feature_cols], function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  })
]

length(numeric_cols)
length(categorical_cols)

numeric_cols <- numeric_cols[!(numeric_cols %in% c("Sex M: 1  F: 2", "ECOG _PS"))]
categorical_cols <- c("Sex M: 1  F: 2", "ECOG _PS", "EGFR_Mutation_Status_mod","IO_Line","Previous_palliative_chemo","Previous_palliative_target","PD-L1_TPS_mod")

get_numeric_cluster_markers <- function(
    analysis_dt,
    cluster_col = "cluster",
    sample_col = "sample",
    numeric_cols,
    min_n_in = 3,
    min_n_out = 3
) {
  res_list <- list()
  clusters <- sort(unique(as.character(analysis_dt[[cluster_col]])))
  
  for (cl in clusters) {
    in_idx <- analysis_dt[[cluster_col]] == cl
    
    for (feat in numeric_cols) {
      x <- analysis_dt[in_idx, get(feat)]
      y <- analysis_dt[!in_idx, get(feat)]
      
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      
      if (length(x) < min_n_in || length(y) < min_n_out) next
      if (length(unique(c(x, y))) < 2) next
      
      wt <- tryCatch(
        wilcox.test(x, y, exact = FALSE),
        error = function(e) NULL
      )
      
      mean_in <- mean(x)
      mean_out <- mean(y)
      var_x <- var(x)
      var_y <- var(y)
      
      sd_pool <- sqrt(
        ((length(x) - 1) * var_x + (length(y) - 1) * var_y) /
          (length(x) + length(y) - 2)
      )
      
      smd <- ifelse(
        is.finite(sd_pool) && sd_pool > 0,
        (mean_in - mean_out) / sd_pool,
        NA_real_
      )
      
      res_list[[length(res_list) + 1]] <- data.table(
        cluster = cl,
        feature = feat,
        n_in = length(x),
        n_out = length(y),
        mean_in = mean_in,
        mean_out = mean_out,
        median_in = median(x),
        median_out = median(y),
        delta_mean = mean_in - mean_out,
        delta_median = median(x) - median(y),
        smd = smd,
        p_value = if (!is.null(wt)) wt$p.value else NA_real_
      )
    }
  }
  
  marker_dt <- rbindlist(res_list, fill = TRUE)
  
  if (nrow(marker_dt) > 0) {
    marker_dt[, fdr := p.adjust(p_value, method = "BH")]
    marker_dt[, direction := fifelse(
      delta_median > 0,
      "high_in_cluster",
      "low_in_cluster"
    )]
    marker_dt[, abs_smd := abs(smd)]
    marker_dt[, score := -log10(fdr + 1e-300) * abs_smd]
    
    setorder(marker_dt, cluster, fdr, -abs_smd)
  }
  
  marker_dt
}

numeric_marker_dt <- get_numeric_cluster_markers(
  analysis_dt = analysis_dt,
  numeric_cols = numeric_cols
)

numeric_marker_dt[, marker_group := fifelse(
  fdr <= 0.1 & smd > 0,
  "high_in_cluster",
  fifelse(
    fdr <= 0.1 & smd < 0,
    "low_in_cluster",
    "not_significant"
  )
)]

ggplot(
  numeric_marker_dt,
  aes(
    x = smd,
    y = -log10(p_value),
    color = marker_group
  )
) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ cluster) +
  scale_color_manual(
    values = c(
      high_in_cluster = "red",
      low_in_cluster = "navy",
      not_significant = "grey70"
    )
  ) +
  theme_classic() +
  labs(
    x = "SMD",
    y = "-log10(p-value)",
    color = "Marker group"
  )#5x9 #20260429_version2_trial1_FAMD_cluster_DEF

ggplot(
  numeric_marker_dt,
  aes(
    x = smd,
    y = -log10(p_value),
    color = marker_group
  )
) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ cluster) +
  scale_color_manual(
    values = c(
      high_in_cluster = "red",
      low_in_cluster = "navy",
      not_significant = "grey70"
    )
  ) +
  geom_text(
    data = numeric_marker_dt[marker_group != "not_significant", ],
    aes(label = feature),
    hjust = -0.1,
    size = 3
  ) +
  theme_classic() +
  labs(
    x = "SMD",
    y = "-log10(p-value)",
    color = "Marker group"
  )#5x9 #20260429_version2_trial1_FAMD_cluster_DEF_v2
get_categorical_cluster_markers <- function(
    analysis_dt,
    cluster_col = "cluster",
    sample_col = "sample",
    categorical_cols,
    min_count = 2
) {
  res_list <- list()
  clusters <- sort(unique(as.character(analysis_dt[[cluster_col]])))
  
  for (cl in clusters) {
    in_idx <- analysis_dt[[cluster_col]] == cl
    
    for (feat in categorical_cols) {
      vals <- as.character(analysis_dt[[feat]])
      vals[is.na(vals)] <- "NA"
      
      levels_i <- sort(unique(vals))
      
      for (lv in levels_i) {
        in_level <- vals == lv
        
        a <- sum(in_idx & in_level)
        b <- sum(in_idx & !in_level)
        c <- sum(!in_idx & in_level)
        d <- sum(!in_idx & !in_level)
        
        if ((a + c) < min_count) next
        
        mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        rownames(mat) <- c("in_cluster", "out_cluster")
        colnames(mat) <- c("level_present", "level_absent")
        
        ft <- tryCatch(
          fisher.test(mat),
          error = function(e) NULL
        )
        
        prop_in <- a / (a + b)
        prop_out <- c / (c + d)
        
        res_list[[length(res_list) + 1]] <- data.table(
          cluster = cl,
          feature = feat,
          level = lv,
          n_in_level = a,
          n_in_total = a + b,
          n_out_level = c,
          n_out_total = c + d,
          prop_in = prop_in,
          prop_out = prop_out,
          delta_prop = prop_in - prop_out,
          odds_ratio = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
          p_value = if (!is.null(ft)) ft$p.value else NA_real_
        )
      }
    }
  }
  
  marker_dt <- rbindlist(res_list, fill = TRUE)
  
  if (nrow(marker_dt) > 0) {
    marker_dt[, fdr := p.adjust(p_value, method = "BH")]
    marker_dt[, direction := fifelse(
      delta_prop > 0,
      "enriched_in_cluster",
      "depleted_in_cluster"
    )]
    marker_dt[, abs_delta_prop := abs(delta_prop)]
    marker_dt[, score := -log10(fdr + 1e-300) * abs_delta_prop]
    
    setorder(marker_dt, cluster, fdr, -abs_delta_prop)
  }
  
  marker_dt
}

categorical_marker_dt <- get_categorical_cluster_markers(
  analysis_dt = analysis_dt,
  categorical_cols = categorical_cols
)

categorical_marker_dt

# 안전하게 copy
cat_plot_dt <- copy(categorical_marker_dt)
# y축 label
cat_plot_dt[, feature_label := paste0(feature, " = ", level)]

# dot color용 그룹
cat_plot_dt[, fdr_group := ifelse(fdr < 0.1, "FDR < 0.1", "NS")]

# 보기 좋게 OR 순으로 정렬
cat_plot_dt <- cat_plot_dt[order(cluster, odds_ratio)]

cat_plot_dt[, feature_label := factor(feature_label, levels = unique(feature_label))]
unique(cat_plot_dt$feature_label)
feature_order <- c('Sex M: 1  F: 2 = 1','Sex M: 1  F: 2 = 2','PD-L1_TPS_mod = NA','PD-L1_TPS_mod = <1','PD-L1_TPS_mod = 1-49','PD-L1_TPS_mod = ≥50','IO_Line = 1st Line','IO_Line = ≥2nd line','EGFR_Mutation_Status_mod = Wild-Type','EGFR_Mutation_Status_mod = Mut','ECOG _PS = 0','ECOG _PS = 1','ECOG _PS = 2','Previous_palliative_chemo = O','Previous_palliative_chemo = X','Previous_palliative_target = O','Previous_palliative_target = X')
cat_plot_dt[, feature_label := factor(feature_label, levels = rev(feature_order))]
ggplot(
  cat_plot_dt,
  aes(
    y = feature_label,
    x = odds_ratio
  )
) +
  geom_segment(
    aes(
      x = 1,
      xend = odds_ratio,
      y = feature_label,
      yend = feature_label
    ),
    color = "grey70"
  ) +
  geom_point(
    aes(color = fdr_group),
    size = 2.5
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_wrap(~ cluster, scales = "free_y") +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "FDR < 0.1" = "red",
      "NS" = "grey50"
    )
  ) +
  theme_classic() +
  labs(
    x = "Odds ratio (log scale)",
    y = "Feature",
    color = ""
  ) #8x16 #20260429_version2_trial1_FAMD_cluster_DEF_categorical
################
# full cohort, single feature add on test
single_scan <- fread('data/20260309_pilot/results/version2/feature_search_base_v2_single_feature_outer/stage2_single_feature_scan.csv')
single_scan <- single_scan[c(order(single_scan$roc_auc_mean)),]
ggplot(single_scan,aes(x=factor(feature_name,levels=single_scan$feature_name), y=roc_auc_mean),group=1)+geom_point(alpha=0.6)+geom_hline(yintercept = 0.6133690476190476,lty='dashed')+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=4))#8x8 #20260429_fullcohort_rocauc
