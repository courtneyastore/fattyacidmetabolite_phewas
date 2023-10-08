library(ggplot2)
library(forcats)
library(tidyverse)
library(dplyr)
library(data.table)
library("ggpubr")
library(plyr)
library(stringr)
library(gridExtra)
library(ggpmisc)

phecode_map_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/phecode_definitions1.2.csv"))
phecode_map_df$phecode <- as.character(phecode_map_df$phecode)
phecode_map_df <- phecode_map_df[c('phecode','phenotype')]

disease_status_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/All_SuggestivePhecodes_O3FA_O6FA_DHA_extractPhecodes.lst"))
disease_status_df$ID <- paste0(disease_status_df$EID,"_",disease_status_df$EID,sep="")
disease_df <- disease_status_df[c('ID','Status_557.1')]
colnames(disease_df) <- c('ID','disease_status')
disease_df <- disease_df[disease_df$disease_status != "Excluded", ]
disease_df$disease_status <- ifelse(disease_df$disease == "Case", "Celiac disease", "Control")

covar_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/12_12_2022_UKBBCovariate_table.tsv"))
covar_df$ID <- paste0(covar_df$eid,"_",covar_df$eid,sep="")
covar_df$sex_corrected <- ifelse(covar_df$sex=="Male",1,0)
covar_df <- covar_df[c('ID','sex','age','age_square','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','sex_corrected')]

environment_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/Zscore_format_normalized_BodyCompositionFieldsUKBB.tsv"))
environment_df$ID <- paste0(environment_df$eid,"_",environment_df$eid,sep="")
environment_df <- environment_df[c('ID','BMI_median','WaistCircumference_median_to_HipCircumference_median_ratio')]

#Remove entries with missing values.
environment_df <- environment_df[!(is.na(environment_df$BMI_median)),]
environment_df <- environment_df[!(is.na(environment_df$WaistCircumference_median_to_HipCircumference_median_ratio)),]

merge_df <- merge(covar_df,disease_df,by="ID")
merge_df <- merge(merge_df,environment_df,by="ID")

merge_df$bmi_group <- ifelse(merge_df$BMI_median > 30,"BMI-Obese","BMI-Non-obese")

female_df <- merge_df[merge_df$sex_corrected == 0, ]
female_df$whr_group <- ifelse(female_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.85,"WHR-Obese","WHR-Non-obese")

male_df <- merge_df[merge_df$sex_corrected == 1, ]
male_df$whr_group <- ifelse(male_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.95,"WHR-Obese","WHR-Non-obese")

final_df <- rbind(male_df,female_df)

colors  <- c("Celiac disease" = "magenta",
             "Control" = "green")

bmi_plt <- ggplot(final_df, aes(x = disease_status, y = BMI_median,fill=disease_status)) +
  geom_boxplot() + scale_fill_manual(values = colors)+
  xlab("Disease status") +
  ylab("BMI") + ggtitle("Celiac disease and BMI\nT-test p-value: 5.74e-39") + 
  theme_bw() + theme(legend.position="none",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))

whr_plt <- ggplot(final_df, aes(x = disease_status, y = WaistCircumference_median_to_HipCircumference_median_ratio,fill=disease_status)) +
  geom_boxplot() +scale_fill_manual(values = colors)+
  xlab("Disease status") +
  ylab("WHR") + ggtitle("Celiac disease and WHR\nT-test p-value: 1.20e-18") + 
  theme_bw() + theme(legend.position="none",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))

ggarrange(bmi_plt, whr_plt, 
          ncol = 2, nrow = 1)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/SuppFig1_CeliacBMIWHR.pdf",height = 8,width = 16,dpi = 300,limitsize = FALSE)

bmi_celiac = t.test(BMI_median ~ disease_status, data = final_df) # 5.739911e-39
whr_celiac = t.test(WaistCircumference_median_to_HipCircumference_median_ratio ~ disease_status, data = final_df)  # 1.201847e-18
