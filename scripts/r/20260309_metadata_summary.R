library(data.table)
library(readxl)
library(ggplot2)
library(reshape2)

meta <- read_excel('data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2.xlsx')
meta <- as.data.frame(meta)
dim(meta[meta$`PD-L1_TPS` != 'NA',])

select.patiant <- fread('data/20260309_pilot/results/20260313_patients_selected.txt', header = F)
meta <- meta[meta$sample_id %in% select.patiant$V1,]

ggplot(meta,aes(fill=`Binarized response`, x=`PFS (Days)`))+geom_histogram()+facet_grid(`Binarized response`~.)
ggplot(meta, aes(x=`Binarized response`,fill=PD_Event,group=PD_Event))+geom_bar()
ggplot(meta, aes(x=`Binarized response`,fill=Death_Event,group=Death_Event))+geom_bar()

meta$`PD-L1_TPS` <- as.numeric(meta$`PD-L1_TPS`) # NA: 16 
meta$`PD-L1_TPS_mod` <- ifelse(is.na(meta$`PD-L1_TPS`),meta$`PD-L1_TPS`,
                               ifelse(meta$`PD-L1_TPS` <1 ,'<1',
                                      ifelse(meta$`PD-L1_TPS` < 50,'1-49','≥50')))
meta$Histology_mod <- ifelse(meta$Histology =='LUAD','LUAD',
                             ifelse(meta$Histology =='LUSC','LUSC','others'))
meta$EGFR_Mutation_Status_mod <- ifelse(meta$EGFR_Mutation_Status =='Wild-Type','Wild-Type','Mut')
# y = Binarized response
meta_melt_binary <- melt(meta, id.vars = c('sample_id','Binarized response'),
                         measure.vars = c('Sex M: 1  F: 2','Age_IO','Histology_mod','Smoking Never: 0 Ex: 1 Current: 2','Drug','ECOG _PS','PFS (Days)','PD_Event','Death_Event','OS (Days)','EGFR_Mutation_Status_mod','IO_Line','Previous_palliative_chemo','Previous_palliative_target','PD-L1_TPS_mod'))
meta_melt_binary$datatype <- ifelse(meta_melt_binary$variable %in% c('Age_IO','PFS (Days)','OS (Days)'),'numeric',
                                    ifelse(meta_melt_binary$variable %in% c('Histology_mod','Smoking Never: 0 Ex: 1 Current: 2','ECOG _PS','Drug','PD-L1_TPS_mod'),'categoric','binary'))
ggplot(meta_melt_binary[meta_melt_binary$datatype == 'categoric',], aes(x=`Binarized response`, fill=value))+geom_bar()+facet_wrap(.~variable)+theme_bw()#4x6 #20260318_catagorical_data
ggplot(meta_melt_binary[meta_melt_binary$datatype == 'binary',], aes(x=`Binarized response`, fill=value))+geom_bar()+facet_wrap(.~variable)+theme_bw()#6x6 #20260318_binary_data

temp <- meta_melt_binary[meta_melt_binary$datatype == 'numeric',]
temp$value <- as.numeric(temp$value)
ggplot(temp[temp$variable != 'Age_IO',], aes(fill=`Binarized response`, x=value))+geom_histogram()+facet_wrap(`Binarized response`~variable)+theme_bw()#4x6 #20260318_numeric_PFS_OS
ggplot(temp[temp$variable == 'Age_IO',], aes(fill=`Binarized response`, x=value))+geom_histogram()+facet_grid(`Binarized response`~.)+theme_bw()#3x4 #20260318_numeric_Age

meta_mod <- meta[,c('sample_id','Binarized response','RECIST response','Study ID','Sample ID','Sex M: 1  F: 2','Age_IO','Histology_mod','Smoking Never: 0 Ex: 1 Current: 2','ECOG _PS','Drug','PD_Event','PFS (Days)','Death_Event','OS (Days)','EGFR_Mutation_Status_mod','IO_Line','Previous_palliative_chemo','Previous_palliative_target','PD-L1_TPS_mod')]
write.table(meta_mod, file='data/20260309_pilot/nsclc_n73/20260309_eQTL Study_SNU (Pilot cohort)-2_mod.txt',row.names = F, col.names = T, sep='\t',quote=F)
# y = PFS
