library(ggplot2)
library(data.table)
library(reshape2)
library(readxl)
library(stringr)

# subtype discovery
v5_importance <- fread('data/20260309_pilot/results/version1/subtype_discovery/base_run_v5/permutation_feature_importance/permutation_importance_summary.csv')
v5_importance.melt <- melt(v5_importance, id.vars=c('feature'), measure.vars = c('base_silhouette','mean_perm_silhouette'))
ggplot(v5_importance.melt, aes(x=factor(feature,levels=v5_importance[c(order(v6_importance$importance)),]$feature), y=value,col=variable))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


v6_importance <- fread('data/20260309_pilot/results/version1/subtype_discovery/base_run_v6/permutation_feature_importance/permutation_importance_summary.csv')
v6_importance[c(order(v6_importance$importance)),]$feature
v6_importance.melt <- melt(v6_importance, id.vars=c('feature'), measure.vars = c('base_silhouette','mean_perm_silhouette'))
ggplot(v6_importance.melt, aes(x=factor(feature,levels=v6_importance[c(order(v6_importance$importance)),]$feature), y=value,col=variable))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


v7_importance <- fread('data/20260309_pilot/results/version1/subtype_discovery/base_run_v7/permutation_feature_importance/permutation_importance_summary.csv')
v7_importance.melt <- melt(v7_importance, id.vars=c('feature'), measure.vars = c('base_silhouette','mean_perm_silhouette'))
ggplot(v7_importance.melt, aes(x=factor(feature,levels=v7_importance[c(order(v7_importance$importance)),]$feature), y=value,col=variable))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ct_importance <- fread('data/20260309_pilot/results/version1/subtype_discovery/pc_only_celltype_sweep/NK_pc1to5/permutation_feature_importance/permutation_importance_summary.csv')
ct_importance.melt <- melt(ct_importance, id.vars=c('feature'), measure.vars = c('base_silhouette_fixed_k','mean_perm_silhouette_fixed_k'))
ggplot(ct_importance.melt, aes(x=factor(feature,levels=ct_importance[c(order(ct_importance$importance_fixed_k)),]$feature), y=value,col=variable))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

sbp_importance <- fread('data/20260309_pilot/results/version1/subtype_discovery/biological_SBP_v1/permutation_feature_importance/permutation_importance_summary.csv')
sbp_importance.melt <- melt(sbp_importance, id.vars=c('feature'), measure.vars = c('base_silhouette_fixed_k','mean_perm_silhouette_fixed_k'))
ggplot(sbp_importance.melt, aes(x=factor(feature,levels=sbp_importance[c(order(sbp_importance$importance_fixed_k)),]$feature), y=value,col=variable))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

v5_pca <- fread('data/20260309_pilot/results/version1/subtype_discovery/base_run_v5/patient_pca.csv')
v5_umap <- fread('data/20260309_pilot/results/version1/subtype_discovery/base_run_v5/patient_umap.csv')
v5_umap$sample <- str_split_fixed(v5_umap$analysis_id,'__',2)[,2]
v5_cluster <- fread('data/20260309_pilot/results/version1/subtype_discovery/base_run_v5/patient_clusters.csv')
v5_umap_cluster <- merge(v5_umap,v5_cluster, by='analysis_id')
meta <- read_excel('data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-1.xlsx')
meta <- as.data.frame(meta)
v5_umap_cluster_meta <- merge(v5_umap_cluster, meta, by.x='sample',by.y='sample_id',all.x=T)
ggplot(v5_umap_cluster_meta,aes(x=UMAP1, y=UMAP2, col=factor(patient_cluster)))+geom_point()
ggplot(v5_umap_cluster_meta,aes(x=UMAP1, y=UMAP2, col=`PFS (Days)`))+geom_point()
ggplot(v5_umap_cluster_meta_pca,aes(x=PC1, y=PC2, col=factor(PD_Event)))+geom_point()
v5_umap_cluster_meta_pca <- merge(v5_umap_cluster_meta, v5_pca, by='analysis_id')
ggplot(v7_umap_cluster_meta_pca,aes(x=PC1, y=PC2, col=factor(patient_cluster)))+geom_point()
ggplot(v7_umap_cluster_meta_pca,aes(x=PC1, y=PC2, col=factor(PD_Event)))+geom_point()

# modeling result
model_metrix <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/model_metrics_summary.csv')
model_metrix <- fread('data/20260309_pilot/results/version1/endpoint_modeling/biological_SBP_v1/model_metrics_summary.csv')
model_metrix$analysis_id <- paste0(model_metrix$analysis_block,'_',model_metrix$endpoint,'_',model_metrix$feature_set)
#roc_auc
#auprc (positive class가 적을때 )
#brier score #예측 확률이 실제 결과와 얼마나 가까운지
#cindex #survival에서 쓰는 순위 성능
ggplot(model_metrix[model_metrix$endpoint=='Binarized_response',], aes(x=feature_set, y=roc_auc_mean))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))#4x5 #20260322_feature_model_binary_sbp_best #20260322_feature_model_binary_sbp #20260318_feature_model_binary
ggplot(model_metrix[model_metrix$endpoint=='PFS',], aes(x=feature_set, y=cindex_mean))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))#4x5 #20260322_feature_model_continuous_sbp_best #20260322_feature_model_continuous_sbp #20260318_feature_model_continuous_sbp
ggplot(model_metrix[model_metrix$endpoint=='PFS_6m_restricted',], aes(x=feature_set, y=cindex_mean))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))#4x5 #20260322_feature_model_continuous_6mRestricted_sbp_best #20260318_feature_model_continuous_sbp
model_metrix.melt <- melt(model_metrix, id.vars=c('endpoint','feature_set'), measure.vars = c('cindex_mean','roc_auc_mean'))
ggplot(model_metrix.melt[!is.na(model_metrix.melt$value),], aes(x=factor(feature_set,levels=c('clinical_only','clinical_plus_bio_sbp','clinical_plus_bio_sbp_plus_pc','clinical_plus_bio_sbp_with_pd_l1','clinical_plus_bio_sbp_plus_pc_with_pd_l1','immune_only_bio_sbp_plus_pc')), y=value))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+facet_grid(.~endpoint)+ylim(c(0.4,0.7))+theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))#4x6 #20260322_feature_model_all_sbp_best
# stability
primary_pfs_clinic <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_only__pfs_stability.csv') #none
primary_binary_clinic <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_only__binarized_response_stability.csv')
primary_binary_immuneAdd <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_plus_prop_plus_pc__binarized_response_stability.csv')
primary_binary_immuneAdd <- primary_binary_immuneAdd[c(order(primary_binary_immuneAdd$selection_frequency)),]
ggplot(primary_binary_immuneAdd, aes(x=factor(feature,levels=primary_binary_immuneAdd$feature), y=selection_frequency))+geom_bar(stat='identity')+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))#4x6 #primary_clinical_plus_prop_plus_pc_binarized_stability
primary_pfs_immuneAdd <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_plus_prop_plus_pc__pfs_stability.csv')
primary_pfs_immuneAdd <- primary_pfs_immuneAdd[c(order(primary_pfs_immuneAdd$selection_frequency)),]
ggplot(primary_pfs_immuneAdd, aes(x=factor(feature,levels=primary_pfs_immuneAdd$feature), y=selection_frequency))+geom_bar(stat='identity')+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))#4x6 #primary_clinical_plus_prop_plus_pc_pfs_stability
secondary_binary_clinic_propAdd <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/secondary__clinical_plus_prop__binarized_response_stability.csv')
ggplot(secondary_binary_clinic_propAdd, aes(x=factor(feature,levels=secondary_binary_clinic_propAdd$feature), y=selection_frequency))+geom_bar(stat='identity')+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))#4x6 #secondary_clinical_plus_prop_binarized_stability
secondary_pfs_clinic_propAdd <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/secondary__clinical_plus_prop__pfs_stability.csv')
secondary_pfs_clinic_propAdd <- secondary_pfs_clinic_propAdd[c(order(secondary_pfs_clinic_propAdd$selection_frequency)),]
ggplot(secondary_pfs_clinic_propAdd, aes(x=factor(feature,levels=secondary_pfs_clinic_propAdd$feature), y=selection_frequency))+geom_bar(stat='identity')+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))#4x6 #secondary_clinical_plus_prop_pfs_stability

# cv_folds
primary_binary_clinic_fold <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_only__binarized_response_cv_folds.csv')
ggplot(primary_binary_clinic_fold, aes(x=alpha,y=roc_auc))+geom_point()
ggplot(primary_binary_clinic_fold, aes(x=l1_ratio,y=roc_auc))+geom_point()
ggplot(primary_binary_clinic_fold, aes(x=inner_auc,y=roc_auc, col=ifelse(alpha ==0.5 & l1_ratio %in% c(0.75,1),'red','black')))+geom_point()+theme_bw() #4x4 #primary_clinical_only_binarized_innerROC_allROC_cor
ggplot(primary_binary_clinic_fold, aes(x=alpha,y=l1_ratio, size=inner_auc))+geom_point(alpha=0.3)+theme_bw() #primary_clinical_only_binarized_hyperparameter_innerROC

primary_pfs_clinic_fold <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_only__pfs_cv_folds.csv')
ggplot(primary_pfs_clinic_fold, aes(x=inner_cindex,y=cindex))+geom_point()+theme_bw() #4x4 #primary_clinical_only_pfs_innerCindex_allCindex_cor
ggplot(primary_pfs_clinic_fold, aes(x=penalizer,y=l1_ratio, size=inner_cindex))+geom_point(alpha=0.3)+theme_bw() #4x5 #primary_clinical_only_pfs_hyperparameter_innerCindex

primary_binary_immuneAdd_fold <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_plus_prop_plus_pc__binarized_response_cv_folds.csv')
ggplot(primary_binary_immuneAdd_fold, aes(x=inner_auc,y=roc_auc))+geom_point()+theme_bw() #4x4 #primary_clinical_plus_prop_plus_pc_binarized_innerROC_allROC_cor
ggplot(primary_binary_immuneAdd_fold, aes(x=alpha,y=l1_ratio, size=inner_auc))+geom_point(alpha=0.3)+theme_bw() #4x5 #primary_clinical_plus_prop_plus_pc_binarized_hyperparameter_innerROC

primary_pfs_immuneAdd_fold <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/primary__clinical_plus_prop_plus_pc__pfs_cv_folds.csv')
ggplot(primary_pfs_immuneAdd_fold, aes(x=inner_cindex,y=cindex))+geom_point()+theme_bw() #4x4 #primary_clinical_plus_prop_plus_pc_pfs_innerCindex_allCindex_cor
ggplot(primary_pfs_immuneAdd_fold, aes(x=penalizer,y=l1_ratio, size=inner_cindex))+geom_point(alpha=0.3)+theme_bw() #4x5 #primary_clinical_plus_prop_plus_pc_pfs_hyperparameter_innerCindex

secondary_binary_prop_fold <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/secondary__clinical_plus_prop__binarized_response_cv_folds.csv')
ggplot(secondary_binary_prop_fold, aes(x=inner_auc,y=roc_auc))+geom_point()+theme_bw() #4x4 #secondary_clinical_plus_prop_binarized_innerROC_allROC_cor
ggplot(secondary_binary_prop_fold, aes(x=alpha,y=l1_ratio, size=inner_auc))+geom_point(alpha=0.3)+theme_bw() #4x5 #secondary_clinical_plus_prop_binarized_hyperparameter_innerROC

secondary_pfs_prop_fold <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/secondary__clinical_plus_prop__pfs_cv_folds.csv')
ggplot(secondary_pfs_prop_fold, aes(x=inner_cindex,y=cindex))+geom_point()+theme_bw() #4x4 #secondary_clinical_plus_prop_pfs_innerCindex_allCindex_cor
ggplot(secondary_pfs_prop_fold, aes(x=penalizer,y=l1_ratio, size=inner_cindex))+geom_point(alpha=0.3)+theme_bw() #4x5 #secondary_clinical_plus_prop_pfs_hyperparameter_innerCindex

# full data fit
full <-fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/full_data_tuning_summary.csv')
full <-fread('data/20260309_pilot/results/version1/endpoint_modeling/biological_SBP_v1/full_data_tuning_summary.csv')
full$analysis_id <- paste0(full$analysis_block,'_',full$endpoint,'_',full$feature_set)
model_metrix_full <- merge(model_metrix, full, by='analysis_id')
model_metrix_full.melt <- melt(model_metrix_full, id.vars=c('analysis_id','endpoint.x','feature_set.x'),measure.vars = c('best_alpha','best_l1_ratio','full_inner_auc','best_penalizer','full_inner_cindex','roc_auc_mean','roc_auc_sd','auprc_mean','auprc_sd','brier_mean','brier_sd','cindex_mean','cindex_sd'))
ggplot(model_metrix_full.melt[model_metrix_full.melt$endpoint.x =='Binarized_response',],
       aes(x=feature_set.x,y=value))+geom_point()+facet_wrap(.~variable, scales='free_y')+theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))#8x8 #binarized_full_data_tuning_summary_sbp_best #binarized_full_data_tuning_summary_sbp #binarized_full_data_tuning_summary
ggplot(model_metrix_full.melt[model_metrix_full.melt$endpoint.x =='PFS',],
       aes(x=feature_set.x,y=value))+geom_point()+facet_wrap(.~variable, scales='free_y')+theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))#8x8 #pfs_full_data_tuning_summary_sbp_best #pfs_full_data_tuning_summary_sbp
ggplot(model_metrix_full.melt[model_metrix_full.melt$endpoint.x =='PFS_6m_restricted',],
       aes(x=feature_set.x,y=value))+geom_point()+facet_wrap(.~variable, scales='free_y')+theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))#8x8 #pfs_6m_full_data_tuning_summary_sbp_best #pfs_full_data_tuning_summary_sbp


full_coeff <- fread('data/20260309_pilot/results/version1/endpoint_modeling/base_run_v4/full_data_coefficients.csv')
full_coeff <- fread('data/20260309_pilot/results/version1/endpoint_modeling/biological_SBP_v1/full_data_coefficients.csv')
ggplot(full_coeff[full_coeff$endpoint=='Binarized_response',], aes(x=feature, y=exp(coef)))+geom_point(alpha=0.5)+facet_wrap(.~feature_set,scale='free_y')+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1,size=5))+geom_hline(yintercept=1,lty='dashed') #6x13 #binarized_coefficient_sbp_best #binarized_coefficient_sbp

ggplot(full_coeff[full_coeff$endpoint=='Binarized_response' & exp(full_coeff$coef) < 25,], aes(x=feature, y=exp(coef)))+geom_point(alpha=0.5)+facet_wrap(.~feature_set,scale='free_y')+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1,size=5))+geom_hline(yintercept=1,lty='dashed')#6x13 #binarized_coefficient_zoom_sbp_best #binarized_coefficient_zoom_sbp

ggplot(full_coeff[full_coeff$endpoint=='PFS',], aes(x=feature, y=exp(coef)))+geom_point(alpha=0.5)+facet_wrap(.~feature_set,scale='free_y')+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1,size=5))+geom_hline(yintercept=1,lty='dashed') #6x13 #pfs_coefficient_sbp_best #pfs_coefficient_sbp

ggplot(full_coeff[full_coeff$endpoint=='PFS_6m_restricted',], aes(x=feature, y=exp(coef)))+geom_point(alpha=0.5)+facet_wrap(.~feature_set,scale='free_y')+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1,size=5))+geom_hline(yintercept=1,lty='dashed') #6x13 #pfs_6mo_coefficient_zoom_best #pfs_coefficient_zoom
