library(ggplot2)
library(forcats)
library(tidyverse)
library(dplyr)
library(data.table)
library("ggpubr")
library(plyr)
library(stringr)
library(gridExtra)


phecode_map_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/phecode_definitions1.2.csv"))
phecode_map_df$phecode <- as.character(phecode_map_df$phecode)
phecode_map_df <- phecode_map_df[c('phecode','phenotype')]

disease_status_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/All_SuggestivePhecodes_O3FA_O6FA_DHA_extractPhecodes.lst"))
disease_lst <- colnames(disease_status_df)[-1]
disease_status_df$ID <- paste0(disease_status_df$EID,"_",disease_status_df$EID,sep="")

prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/Zscore_format_normalized_metabolite_table2_participant.tsv"))
prs_df$ID <- paste0(prs_df$eid,"_",prs_df$eid,sep="")
prs_df <- prs_df[c('ID','p23444_zscore')]
prs_df <- prs_df[!is.na(prs_df$p23444_zscore), ]
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
#merge_df <- merge(merge_df,disease_status_df,by="ID")
merge_df <- merge(merge_df,environment_df,by="ID")
merge_df$bmi_group <- ifelse(merge_df$BMI_median > 30,1,0)

female_df <- merge_df[merge_df$sex_corrected == 0, ]
female_df$whr_group <- ifelse(female_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.85,1,0)

male_df <- merge_df[merge_df$sex_corrected == 1, ]
male_df$whr_group <- ifelse(male_df$WaistCircumference_median_to_HipCircumference_median_ratio > 0.95,1,0)

merge_df <- rbind(male_df,female_df)

#disease_lst <- disease_lst[1]
model_colors <- c("BMI-Obese" = "#1f9990","BMI-Non-obese" = "#CCFFFF","WHR-Obese"="#9F2B68","WHR-Non-obese"="#F8C8DC")

for (disease in disease_lst){
  
  phecode = gsub("Status_","",disease)
  
  disease_phecode_df <- phecode_map_df[phecode_map_df$phecode == phecode, ]
  
  disease_name = disease_phecode_df$phenotype[1]
  
  #disease_name_file = gsub(" ","_",disease_name)
  disease_name_file = gsub("[^A-Za-z0-9]", "_", disease_name)
  output_file_name = paste0("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/LevelsO3FA_BMI_WHR_PrevPerc_Batch/LevelsO3FA_BMI_WHR_",disease_name_file,".pdf")
  print(output_file_name)
  
  disease_name_prev_label = paste0("Prevalence of ",disease_name," (%)")
  
  disease_df <- disease_status_df[c('ID',disease)]
  colnames(disease_df) <- c('ID','disease_status')
  
  disease_df <- disease_df[disease_df$disease_status != "Excluded", ]

  disease_merge_disease_df <- merge(merge_df,disease_df,by="ID")
  disease_merge_disease_df$disease_status <- ifelse(disease_merge_disease_df$disease == "Case", 1, 0)
  
  # Obese BMI
  obese_bmi_df <- disease_merge_disease_df[disease_merge_disease_df$bmi_group == 1, ]
  obese_bmi_df <- obese_bmi_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
  model1_prev_perc_df <- obese_bmi_df %>%
    group_by(Percentile_PRS) %>%
    dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
  model1_prev_perc_df <- as.data.frame(model1_prev_perc_df)
  model1_prev_perc_df$model <- "BMI-Obese"
  colnames(model1_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
  model1_prev_perc_df <- model1_prev_perc_df[c('PRS_percentile','Prev','model')]
  
  # Non-obese BMI
  non_obese_bmi_df <- disease_merge_disease_df[disease_merge_disease_df$bmi_group == 0, ]
  non_obese_bmi_df <- non_obese_bmi_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
  model2_prev_perc_df <- non_obese_bmi_df %>%
    group_by(Percentile_PRS) %>%
    dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
  model2_prev_perc_df <- as.data.frame(model2_prev_perc_df)
  model2_prev_perc_df$model <- "BMI-Non-obese"
  colnames(model2_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
  model2_prev_perc_df <- model2_prev_perc_df[c('PRS_percentile','Prev','model')]
  
  # Obese WHR
  obese_whr_df <- disease_merge_disease_df[disease_merge_disease_df$whr_group == 1, ]
  obese_whr_df <- obese_whr_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
  model3_prev_perc_df <- obese_whr_df %>%
    group_by(Percentile_PRS) %>%
    dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
  model3_prev_perc_df <- as.data.frame(model3_prev_perc_df)
  model3_prev_perc_df$model <- "WHR-Obese"
  colnames(model3_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
  model3_prev_perc_df <- model3_prev_perc_df[c('PRS_percentile','Prev','model')]
  
  # Non-obese WHR
  non_obese_whr_df <- disease_merge_disease_df[disease_merge_disease_df$whr_group == 0, ]
  non_obese_whr_df <- non_obese_whr_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
  model4_prev_perc_df <- non_obese_whr_df %>%
    group_by(Percentile_PRS) %>%
    dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
  model4_prev_perc_df <- as.data.frame(model4_prev_perc_df)
  model4_prev_perc_df$model <- "WHR-Non-obese"
  colnames(model4_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
  model4_prev_perc_df <- model4_prev_perc_df[c('PRS_percentile','Prev','model')]
  
  final_df <- rbind(model1_prev_perc_df,model2_prev_perc_df,model3_prev_perc_df,model4_prev_perc_df)
  
  plt <- ggplot(final_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
    theme_bw() + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    xlab("Percentile of Omega-3 fatty acid levels") + ylab(disease_name_prev_label) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) +
    scale_fill_manual(name="Group",values=model_colors)+ 
    scale_color_manual(name="Group",values=model_colors)
  
  ggsave(output_file_name,plot = plt,height = 12,width = 12,dpi = 300,limitsize = FALSE)
  rm(plt)
}
