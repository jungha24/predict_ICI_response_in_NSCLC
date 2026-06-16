## date: 2026/04/03
## version2 base

## library
library(data.table)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(ggbeeswarm)

## path
model_outdir = 'data/20260309_pilot/results/version2/feature_search_base'
model_outdir_v2 = 'data/20260309_pilot/results/version2/feature_search_base_v2'
## clincial only
file_path <- file.path(model_outdir, "stage2_baseline_metrics.csv")
baseline_mat <- fread(file_path)

## immune feature add on file
## 1) outer fold per row
## 1-1) binarized response
file_path <- file.path(model_outdir, "stage2_cv_details__binarized_response.csv")
stage2_cv_detail_binary <- fread(file_path)
df_mean <- stage2_cv_detail_binary[stage2_cv_detail_binary$is_baseline ==FALSE,] %>%
  group_by(feature_name) %>%
  summarise(mean_auc = mean(roc_auc, na.rm = TRUE))
head(df_mean[c(order(-df_mean$mean_auc)),])
df_mean <- df_mean[c(order(-df_mean$mean_auc)),]

ggplot(df_mean, aes(x=factor(feature_name,levels=df_mean$feature_name), y=mean_auc))+
  geom_point(size=0.8)+
  geom_hline(yintercept = mean(stage2_cv_detail_binary[stage2_cv_detail_binary$is_baseline==TRUE,]$roc_auc),lty='dashed',col='pink')+
  theme_few()+theme(axis.text.x = element_text(size=1, hjust=1, vjust=1,angle=90))#4x9 #20260403_stage2_cv_detail__binarized_response_single_feature_performance
## 1-2) PFS
file_path <- file.path(model_outdir, "stage2_cv_details__pfs.csv")
stage2_cv_detail_pfs <- fread(file_path)
df_mean <- stage2_cv_detail_pfs[stage2_cv_detail_pfs$is_baseline ==FALSE,] %>%
  group_by(feature_name) %>%
  summarise(mean_cindex = mean(cindex, na.rm = TRUE))
head(df_mean[c(order(-df_mean$mean_cindex)),])
df_mean <- df_mean[c(order(-df_mean$mean_cindex)),]
ggplot(df_mean, aes(x=factor(feature_name,levels=df_mean$feature_name), y=mean_cindex))+
  geom_point(size=0.8)+
  geom_hline(yintercept = mean(stage2_cv_detail_pfs[stage2_cv_detail_pfs$is_baseline==TRUE,]$cindex),lty='dashed',col='pink')+
  theme_few()+theme(axis.text.x = element_text(size=1, hjust=1, vjust=1,angle=90))#4x9 #20260403_stage2_cv_detail__pfs_single_feature_performance
## 1-3) PFS 6mo restricted
file_path <- file.path(model_outdir, "stage2_cv_details__pfs_6m_restricted.csv")
stage2_cv_detail_pfs_6m <- fread(file_path)
df_mean <- stage2_cv_detail_pfs_6m[stage2_cv_detail_pfs_6m$is_baseline ==FALSE,] %>%
  group_by(feature_name) %>%
  summarise(mean_cindex = mean(cindex, na.rm = TRUE))
head(df_mean[c(order(-df_mean$mean_cindex)),])
df_mean <- df_mean[c(order(-df_mean$mean_cindex)),]
ggplot(df_mean, aes(x=factor(feature_name,levels=df_mean$feature_name), y=mean_cindex))+
  geom_point(size=0.8)+
  geom_hline(yintercept = mean(stage2_cv_detail_pfs_6m[stage2_cv_detail_pfs_6m$is_baseline==TRUE,]$cindex),lty='dashed',col='pink')+
  theme_few()+theme(axis.text.x = element_text(size=1, hjust=1, vjust=1,angle=90))#4x9 #20260403_stage2_cv_detail__pfs_6mo_restricted_single_feature_performance
## 2) event per test ratio per fold
event_size_pfs_6mo <- stage2_cv_detail_pfs_6m %>%
  group_by(fold) %>%
  summarise(
    event_rate = sum(events_test, na.rm = TRUE) / sum(n_test, na.rm = TRUE)
  )
event_size_pfs_6mo$endpoint <- 'pfs_6mo_restricted'
event_size_pfs <- stage2_cv_detail_pfs %>%
  group_by(fold) %>%
  summarise(
    event_rate = sum(events_test, na.rm = TRUE) / sum(n_test, na.rm = TRUE)
  )
event_size_pfs$endpoint <- 'pfs'
event_size_binary <- stage2_cv_detail_binary %>%
  group_by(fold) %>%
  summarise(
    event_rate = sum(events_test, na.rm = TRUE) / sum(n_test, na.rm = TRUE)
  )
event_size_binary$endpoint <- 'binarized_response'
event_size <- rbind(event_size_binary, event_size_pfs, event_size_pfs_6mo)
ggplot(event_size, aes(x=endpoint, y= event_rate))+geom_boxplot(outlier.size = 0)+geom_beeswarm(alpha=0.6,size=1)+theme_bw()+ylim(c(0,1))#4x4 #20260403_event_n_per_test_n_per_fold
## 3) summarize
file_path <- file.path(model_outdir, "stage2_single_feature_scan.csv")
stage2_scan <- fread(file_path)
ggplot(stage2_scan[stage2_scan$endpoint_name=='Binarized_response',], aes(x=feature_name, y=roc_auc_mean))+geom_point()
# temp <- stage2_scan[stage2_scan$endpoint_name=='Binarized_response',]
temp <- stage2_scan[stage2_scan$endpoint_name=='PFS_6m_restricted',]
# temp <- temp[c(order(-temp$roc_auc_mean))]
temp <- temp[c(order(-temp$cindex_mean))]
# ggplot(temp[temp$delta_roc_auc_mean >0,], aes(x = factor(feature_name,levels=temp$feature_name), y = roc_auc_mean)) +
#   geom_line(aes(group = 1)) +
#   geom_point() +
#   geom_errorbar(aes(
#     ymin = roc_auc_mean - roc_auc_sd,
#     ymax = roc_auc_mean + roc_auc_sd
#   ),
#   width = 0.2
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6)
#   ) +
#   labs(
#     x = "feature_name",
#     y = "roc_auc_mean"
#   )#6x10 #20260403_stage2_binarized_response_delta_positive
# ggplot(temp[temp$delta_roc_auc_mean >0,], aes(x = factor(feature_name,levels=temp$feature_name), y = full_data_coef)) +
#   geom_line(aes(group = 1)) +
#   geom_point() +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6)
#   ) +
#   labs(
#     x = "feature_name",
#     y = "full_data_coef"
#   )+geom_hline(yintercept = 0,lty='dashed')#6x6 #20260403_stage2_binarized_response_delta_positive_fulldata_coeff
# ggplot(temp[temp$delta_cindex_mean >0,], aes(x = factor(feature_name,levels=temp$feature_name), y = cindex_mean)) +
#   geom_line(aes(group = 1)) +
#   geom_point() +
#   geom_errorbar(aes(
#     ymin = cindex_mean - cindex_sd,
#     ymax = cindex_mean + cindex_sd
#   ),
#   width = 0.2
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6)
#   ) +
#   labs(
#     x = "feature_name",
#     y = "cindex_mean"
#   )#6x10 #20260403_stage2_pfs_6mo_delta_positive
ggplot(temp[temp$delta_cindex_mean >0,], aes(x = factor(feature_name,levels=temp$feature_name), y = full_data_coef)) +
  geom_line(aes(group = 1)) +
  geom_point(size=1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6)
  ) +
  labs(
    x = "feature_name",
    y = "full_data_coef"
  )+geom_hline(yintercept = 0,lty='dashed')#6x10 #20260403_stage2_pfs_6mo_delta_positive_fulldata_coeff

### feature pruning
file_path <- file.path(model_outdir, "stage3_pruned_candidates.csv")
pruned_cand <- fread(file_path)
file_path <- file.path(model_outdir_v2, "stage3_pruned_candidates.csv")
pruned_cand_v2 <- fread(file_path)

## stage4
## 1-1) binarized response
file_path <- file.path(model_outdir, "stage4_top10__binary__Binarized_response.csv")
stage4_binary_subset <- fread(file_path)
## 1-2) PFS
file_path <- file.path(model_outdir, "stage4_top10__survival__PFS.csv")
stage4_binary_subset <- fread(file_path)
## 1-3) PFS_6mo
file_path <- file.path(model_outdir, "stage4_top10__survival__PFS_6m_restricted.csv")
stage4_binary_subset <- fread(file_path)
## family cap
file_path <- file.path(model_outdir, "feature_catalog_resolved.csv")
family <- fread(file_path)

## stage 5
file_path <- file.path(model_outdir, "stage5_best_subsets.csv")
result <- fread(file_path)
file_path <- file.path(model_outdir, "stage5_best_subset_bootstrap_stability.csv")
result_bs <- fread(file_path)




## date: 2026/04/06
## version2 base
## version2-v2 + base_stage2 부터

## library
library(data.table)
library(ggplot2)
library(ggthemes)
## path
model_outdir = 'data/20260309_pilot/results/version2/feature_search_base'
model_outdir_v2 = 'data/20260309_pilot/results/version2/feature_search_base_v2'

## clincial only
file_path <- file.path(model_outdir, "stage2_baseline_metrics.csv")
baseline_mat <- fread(file_path)
clinical_only <- baseline_mat[baseline_mat$endpoint_name=='Binarized_response',]$roc_auc_mean

## stage4 best subset
file_path <- file.path("data/20260309_pilot/results/version2/feature_search_base_v2/start_from_base_v1_stage2/stage4_subset_results.csv")
stage4_all <- fread(file_path)
stage4_all$feature <- paste0('feature_comb',seq(1,nrow(stage4_all)))
stage4_all <- stage4_all[c(order(-stage4_all$roc_auc_mean)),]
ggplot(stage4_all, aes(x=factor(feature,levels=stage4_all$feature), y=roc_auc_mean))+geom_point(size=0.6)+geom_line()+geom_hline(yintercept = clinical_only, lty='dashed',col='pink')+theme_few()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=2))#4x8 #20260406_v2_start_from_stage2_stage4_all
file_path <- file.path("data/20260309_pilot/results/version2/feature_search_base_v2/start_from_base_v1_stage2/stage4_top10__binary__Binarized_response.csv")
stage4_best <- fread(file_path)
stage4_best$feature <- paste0('feature_comb',seq(1,nrow(stage4_best)))
stage4_best <- stage4_best[c(order(-stage4_best$roc_auc_mean)),]
ggplot(stage4_best, aes(x=factor(feature,levels=stage4_best$feature), y=roc_auc_mean))+geom_bar(stat='identity',fill='grey60')+geom_hline(yintercept = clinical_only, lty='dashed',col='pink')+theme_few()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))#4x3 (p) #20260406_v2_start_from_stage2_stage4_top10

## version2 
## stage2
file_path <- file.path("data/20260309_pilot/results/version2/feature_search_base_v2/stage2_baseline_metrics.csv")
baseline_mat <- fread(file_path)
with_pd_l1 <- baseline_mat[baseline_mat$endpoint_name=='Binarized_response',]$roc_auc_mean

file_path <- file.path("data/20260309_pilot/results/version2/feature_search_base_v2/stage2_single_feature_scan.csv")
stage2_all <- fread(file_path)
stage2_all <- stage2_all[c(order(-stage2_all$roc_auc_mean)),]
ggplot(stage2_all, aes(x=factor(feature_name,levels=stage2_all$feature_name), y=roc_auc_mean))+geom_point(size=0.6)+geom_line()+geom_hline(yintercept = with_pd_l1, lty='dashed',col='pink')+theme_few()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=2))#4x9 #20260406_version2_with_pdl1tps_stage2_cv_detail__binarized_response_single_feature_performance

file_path <- file.path("data/20260309_pilot/results/version2/feature_search_base_v2/stage4_top10__binary__Binarized_response.csv")
stage4_best <- fread(file_path)
stage4_best$feature <- paste0('feature_comb',seq(1,nrow(stage4_best)))
stage4_best <- stage4_best[c(order(-stage4_best$roc_auc_mean)),]
ggplot(stage4_best, aes(x=factor(feature,levels=stage4_best$feature), y=roc_auc_mean))+geom_bar(stat='identity',fill='grey60')+geom_hline(yintercept = with_pd_l1, lty='dashed',col='pink')+theme_few()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))#4x3 (p) #20260406_version2_with_pdl1tps_stage4_top10

library('readxl')
year <- c(2025,2024,2023,2022)
n <- c(68,62,92,115)
temp <- data.frame(year,n)
ggplot(temp,aes(x=factor(year,levels=c(2022,2023,2024,2025)),y=n))+geom_bar(stat='identity')+theme_few()#3x3 #20260406_n_clinical_ICI
