library(ggplot2)
library(forcats)
library(tidyverse)
library(dplyr)
library(data.table)
library("ggpubr")
library(plyr)
library(stringr)
library(gridExtra)
library(ggrepel)
library(rsq)
library(matrixStats)
options(ggrepel.max.overlaps = Inf)


disease_status_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/V2ALLmerged_extractPhecodes.lst"))
disease_status_df$ID <- paste0(disease_status_df$EID,"_",disease_status_df$EID,sep="")
disease_status_df <- disease_status_df[c('ID','Status_550.2')]

prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmergedPRS_UKBB.tsv"))
prs_df <- prs_df[c('ID','PRS_met-d-omega6_PRS')]
colnames(prs_df) <- c('ID','PRS')

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

merge_df <- merge(covar_df,prs_df,by="ID")
merge_df <- merge(merge_df,disease_status_df,by="ID")
merge_df <- merge(merge_df,environment_df,by="ID")
merge_df <- merge_df[merge_df$Status_550.2 != "Excluded", ]
merge_df$disease_status <- ifelse(merge_df$Status_550.2 == "Case", 1, 0)
merge_df$bmi_group <- ifelse(merge_df$BMI_median > 30,1,0)

female_df <- merge_df[merge_df$sex_corrected == 0, ]
female_df$whr_group <- ifelse(female_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.85,1,0)

male_df <- merge_df[merge_df$sex_corrected == 1, ]
male_df$whr_group <- ifelse(male_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.95,1,0)

merge_df <- rbind(male_df,female_df)

# Obese BMI
obese_bmi_df <- merge_df[merge_df$bmi_group == 1, ]
obese_bmi_df <- obese_bmi_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model1_prev_perc_df <- obese_bmi_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model1_prev_perc_df <- as.data.frame(model1_prev_perc_df)
model1_prev_perc_df$model <- "BMI-Obese"
colnames(model1_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model1_prev_perc_df <- model1_prev_perc_df[c('PRS_percentile','Prev','model')]

# Non-obese BMI
non_obese_bmi_df <- merge_df[merge_df$bmi_group == 0, ]
non_obese_bmi_df <- non_obese_bmi_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model2_prev_perc_df <- non_obese_bmi_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model2_prev_perc_df <- as.data.frame(model2_prev_perc_df)
model2_prev_perc_df$model <- "BMI-Non-obese"
colnames(model2_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model2_prev_perc_df <- model2_prev_perc_df[c('PRS_percentile','Prev','model')]

# Obese WHR
obese_whr_df <- merge_df[merge_df$whr_group == 1, ]
obese_whr_df <- obese_whr_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model3_prev_perc_df <- obese_whr_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model3_prev_perc_df <- as.data.frame(model3_prev_perc_df)
model3_prev_perc_df$model <- "WHR-Obese"
colnames(model3_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model3_prev_perc_df <- model3_prev_perc_df[c('PRS_percentile','Prev','model')]

# Non-obese WHR
non_obese_whr_df <- merge_df[merge_df$whr_group == 0, ]
non_obese_whr_df <- non_obese_whr_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model4_prev_perc_df <- non_obese_whr_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model4_prev_perc_df <- as.data.frame(model4_prev_perc_df)
model4_prev_perc_df$model <- "WHR-Non-obese"
colnames(model4_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model4_prev_perc_df <- model4_prev_perc_df[c('PRS_percentile','Prev','model')]


final_df <- rbind(model1_prev_perc_df,model2_prev_perc_df,model3_prev_perc_df,model4_prev_perc_df)

model_colors <- c("BMI-Obese" = "#1f9990","BMI-Non-obese" = "#CCFFFF","WHR-Obese"="#9F2B68","WHR-Non-obese"="#F8C8DC")

ggplot(final_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of PGS-Omega-6 fatty acids") + ylab("Prevalence of Diaphragmatic hernia (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) +
  scale_fill_manual(name="Group",values=model_colors)+ 
  scale_color_manual(name="Group",values=model_colors)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Figure6/V2_O6FA_DiaphragmaticHernia_BMI_WHR_PGSPrev_plt.pdf",height = 8,width = 8,dpi = 300,limitsize = FALSE)













