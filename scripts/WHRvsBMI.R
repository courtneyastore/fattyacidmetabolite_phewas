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

merge_df <- merge(covar_df,environment_df,by="ID")

merge_df <- merge_df[merge_df$WaistCircumference_median_to_HipCircumference_median_ratio < 2, ]
merge_df$bmi_group <- ifelse(merge_df$BMI_median > 30,"BMI-Obese","BMI-Non-obese")

female_df <- merge_df[merge_df$sex_corrected == 0, ]
female_df$whr_group <- ifelse(female_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.85,"WHR-Obese","WHR-Non-obese")

male_df <- merge_df[merge_df$sex_corrected == 1, ]
male_df$whr_group <- ifelse(male_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.95,"WHR-Obese","WHR-Non-obese")

final_df <- rbind(male_df,female_df)

final_df$obesity_status <- ifelse(final_df$bmi_group == "BMI-Obese" & final_df$whr_group == "WHR-Obese", "WHR-Obese and BMI-Obese",
                                  ifelse(final_df$bmi_group == "BMI-Non-obese" & final_df$whr_group == "WHR-Obese", "WHR-Obese and BMI-Non-obese",
                                         ifelse(final_df$bmi_group == "BMI-Obese" & final_df$whr_group == "WHR-Non-obese", "WHR-Non-obese and BMI-Obese",
                                                ifelse(final_df$bmi_group == "BMI-Non-obese" & final_df$whr_group == "WHR-Non-obese", "WHR-Non-obese and BMI-Non-obese", NA))))

colors  <- c("WHR-Non-obese and BMI-Non-obese" = "lightblue",
             "WHR-Non-obese and BMI-Obese" = "purple",
             "WHR-Obese and BMI-Non-obese" = "blue",
             "WHR-Obese and BMI-Obese"="red")

ggplot(final_df, aes(x=BMI_median, y=WaistCircumference_median_to_HipCircumference_median_ratio,color=obesity_status)) +
  scale_color_manual(values = colors)+ 
  geom_point(size=2,alpha=0.1) + ylab("WHR") + xlab("BMI") + 
  theme_bw() + theme(legend.position="none",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/SuppFig1_WHR_vs_BMI_obesity_thresholds.jpg",height = 8,width = 8,dpi = 300,limitsize = FALSE)

ggplot(final_df, aes(x=BMI_median, y=WaistCircumference_median_to_HipCircumference_median_ratio)) +
  geom_point(size=2,alpha=0.1) + ylab("WHR") + xlab("BMI") + 
  theme_bw() + theme(legend.position="none",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/SuppFig1_WHR_vs_BMI.jpg",height = 8,width = 8,dpi = 300,limitsize = FALSE)

ggplot(final_df, aes(x=BMI_median, y=WaistCircumference_median_to_HipCircumference_median_ratio)) +
  ylab("WHR") + xlab("BMI") + geom_bin_2d(bins=200) + 
  scale_fill_continuous(type = "viridis") +
  theme_bw() + theme(legend.key.width = unit(1, 'in'),legend.key.height= unit(0.5, 'in'),legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/SuppFig1_WHR_vs_BMI_density.jpg",height = 8,width = 8,dpi = 300,limitsize = FALSE)











bmi_celiac = t.test(BMI_median ~ disease_status, data = final_df) # 5.739911e-39
whr_celiac = t.test(WaistCircumference_median_to_HipCircumference_median_ratio ~ disease_status, data = final_df)  # 1.201847e-18
