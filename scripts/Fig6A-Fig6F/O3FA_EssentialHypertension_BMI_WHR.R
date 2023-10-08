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
disease_status_df <- disease_status_df[c('ID','Status_401.1')]

prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmergedPRS_UKBB.tsv"))
prs_df <- prs_df[c('ID','PRS_met-d-omega3_PRS')]
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
merge_df <- merge_df[merge_df$Status_401.1 != "Excluded", ]
merge_df$disease_status <- ifelse(merge_df$Status_401.1 == "Case", 1, 0)
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
  xlab("Percentile of PGS-Omega-3 fatty acids") + ylab("Prevalence of Essential hypertension (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) +
  scale_fill_manual(name="Group",values=model_colors)+ 
  scale_color_manual(name="Group",values=model_colors)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Figure6/V2_O3FA_EssentialHypertension_BMI_WHR_PGSPrev_plt.pdf",height = 8,width = 8,dpi = 300,limitsize = FALSE)













metabolite_environment_df <- merge(metabolite_df,environment_df,by="ID")
MASTER_metabolite_environment_df <- merge(merge_df,metabolite_environment_df,by="ID")

# ------------------- Assess high vs. low DHA for Depression -----------------------
MASTER_metabolite_df$p23444_zscore_group <- ifelse(MASTER_metabolite_df$p23444_zscore >= mean(MASTER_metabolite_df$p23444_zscore),1,0) # 1 is HIGH DHA and 0 id LOW DHA

metabolite_gxe_glm = glm(disease_status~PRS_ieu_a_1188*p23444_zscore_group+PRS_ieu_a_1188+p23444_zscore_group+sex_corrected+age+age_square+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=MASTER_metabolite_df,family="binomial")
summary(metabolite_gxe_glm)

# Get low DHA
low_dha_df <- MASTER_metabolite_df[MASTER_metabolite_df$p23444_zscore_group == 0, ]
low_dha_df <- low_dha_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model1_prev_perc_df <- low_dha_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model1_prev_perc_df <- as.data.frame(model1_prev_perc_df)
model1_prev_perc_df$model <- "Low-Omega-3 fatty acids"
colnames(model1_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model1_prev_perc_df <- model1_prev_perc_df[c('PRS_percentile','Prev','model')]

# Get high DHA
high_dha_df <- MASTER_metabolite_df[MASTER_metabolite_df$p23444_zscore_group == 1, ]
high_dha_df <- high_dha_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model2_prev_perc_df <- high_dha_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model2_prev_perc_df <- as.data.frame(model2_prev_perc_df)
model2_prev_perc_df$model <- "High-Omega-3 fatty acids"
colnames(model2_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model2_prev_perc_df <- model2_prev_perc_df[c('PRS_percentile','Prev','model')]

low_high_DHA_prs_prev_df <- rbind(model1_prev_perc_df,model2_prev_perc_df)

model_colors <- c("Low-Omega-3 fatty acids" = "royalblue1","High-Omega-3 fatty acids" = "blue4")

low_high_o3fa_prev_plt <- ggplot(low_high_DHA_prs_prev_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of PGS-Major depressive disorder") + ylab("Prevalence of disease (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) +
  scale_fill_manual(name="Group",values=model_colors)+ ylim(0,20) +
  scale_color_manual(name="Group",values=model_colors)
low_high_o3fa_prev_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Figure6/O3FA_MajorDepression_PGSPrev_plt.pdf",height = 8,width = 8,dpi = 300,limitsize = FALSE)


# ------------------- Assess high vs. low BMI for Depression -----------------------
MASTER_environment_df$Zscore_BMI_median_group <- ifelse(MASTER_environment_df$Zscore_BMI_median >= mean(MASTER_environment_df$Zscore_BMI_median),1,0) # 1 is HIGH DHA and 0 id LOW DHA

environment_gxe_glm = glm(disease_status~PRS_ieu_a_1188*Zscore_BMI_median_group+PRS_ieu_a_1188+Zscore_BMI_median_group+sex_corrected+age+age_square+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=MASTER_environment_df,family="binomial")
summary(environment_gxe_glm)

# Get low BMI
low_bmi_df <- MASTER_environment_df[MASTER_environment_df$Zscore_BMI_median_group == 0, ]
low_bmi_df <- low_bmi_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model3_prev_perc_df <- low_bmi_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model3_prev_perc_df <- as.data.frame(model3_prev_perc_df)
model3_prev_perc_df$model <- "Low-BMI"
colnames(model3_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model3_prev_perc_df <- model3_prev_perc_df[c('PRS_percentile','Prev','model')]

# Get high BMI
high_bmi_df <- MASTER_environment_df[MASTER_environment_df$Zscore_BMI_median_group == 1, ]
high_bmi_df <- high_bmi_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model4_prev_perc_df <- high_bmi_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model4_prev_perc_df <- as.data.frame(model4_prev_perc_df)
model4_prev_perc_df$model <- "High-BMI"
colnames(model4_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model4_prev_perc_df <- model4_prev_perc_df[c('PRS_percentile','Prev','model')]

low_high_BMI_prs_prev_df <- rbind(model3_prev_perc_df,model4_prev_perc_df)

model_colors <- c("High-BMI" = "#6C0BA9","Low-BMI" = "#D7A1F9")

low_high_bmi_prev_plt <- ggplot(low_high_BMI_prs_prev_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of PGS-Major depressive disorder") + ylab("Prevalence of disease (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) + ylim(0,20) +
  scale_fill_manual(name="Group",values=model_colors)+
  scale_color_manual(name="Group",values=model_colors)
low_high_bmi_prev_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Figure6/O3FA_BMI_PGSPrev_plt.pdf",height = 8,width = 8,dpi = 300,limitsize = FALSE)

# ------------------- Assess high vs. low DHA and BMI for Depression -----------------------
MASTER_metabolite_environment_df$p23444_zscore_group <- ifelse(MASTER_metabolite_environment_df$p23444_zscore >= mean(MASTER_metabolite_environment_df$p23444_zscore),1,0) # 1 is HIGH DHA and 0 id LOW DHA
MASTER_metabolite_environment_df$Zscore_BMI_median_group <- ifelse(MASTER_metabolite_environment_df$Zscore_BMI_median >= mean(MASTER_metabolite_environment_df$Zscore_BMI_median),1,0) # 1 is HIGH DHA and 0 id LOW DHA

# LOW LOW
lowENV1_lowENV2_df <- MASTER_metabolite_environment_df[(MASTER_metabolite_environment_df$Zscore_BMI_median_group==0) & (MASTER_metabolite_environment_df$p23444_zscore_group==0),]
lowENV1_lowENV2_df <- lowENV1_lowENV2_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model5_prev_perc_df <- lowENV1_lowENV2_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model5_prev_perc_df <- as.data.frame(model5_prev_perc_df)
model5_prev_perc_df$model <- "Low-Omega-3 fatty acids & Low-BMI"
colnames(model5_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model5_prev_perc_df <- model5_prev_perc_df[c('PRS_percentile','Prev','model')]

# HIGH HIGH
highENV1_highENV2_df <- MASTER_metabolite_environment_df[(MASTER_metabolite_environment_df$Zscore_BMI_median_group==1) & (MASTER_metabolite_environment_df$p23444_zscore_group==1),]
highENV1_highENV2_df <- highENV1_highENV2_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model6_prev_perc_df <- highENV1_highENV2_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model6_prev_perc_df <- as.data.frame(model6_prev_perc_df)
model6_prev_perc_df$model <- "High-Omega-3 fatty acids & High-BMI"
colnames(model6_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model6_prev_perc_df <- model6_prev_perc_df[c('PRS_percentile','Prev','model')]

# LOW HIGH
lowENV1_highENV2_df <- MASTER_metabolite_environment_df[(MASTER_metabolite_environment_df$Zscore_BMI_median_group==0) & (MASTER_metabolite_environment_df$p23444_zscore_group==1),]
lowENV1_highENV2_df <- lowENV1_highENV2_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model7_prev_perc_df <- lowENV1_highENV2_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model7_prev_perc_df <- as.data.frame(model7_prev_perc_df)
model7_prev_perc_df$model <- "High-Omega-3 fatty acids & Low-BMI"
colnames(model7_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model7_prev_perc_df <- model7_prev_perc_df[c('PRS_percentile','Prev','model')]

# HIGH LOW
highENV1_lowENV2_df <- MASTER_metabolite_environment_df[(MASTER_metabolite_environment_df$Zscore_BMI_median_group==1) & (MASTER_metabolite_environment_df$p23444_zscore_group==0),]
highENV1_lowENV2_df <- highENV1_lowENV2_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model8_prev_perc_df <- highENV1_lowENV2_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model8_prev_perc_df <- as.data.frame(model8_prev_perc_df)
model8_prev_perc_df$model <- "Low-Omega-3 fatty acids & High-BMI"
colnames(model8_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model8_prev_perc_df <- model8_prev_perc_df[c('PRS_percentile','Prev','model')]

env1_env2_prs_prev_df <- rbind(model5_prev_perc_df,model6_prev_perc_df,model7_prev_perc_df,model8_prev_perc_df)

model_colors <- c("Low-Omega-3 fatty acids & High-BMI" = "#DB4437","High-Omega-3 fatty acids & Low-BMI" = "#4285F4","High-Omega-3 fatty acids & High-BMI"="#F4B400","Low-Omega-3 fatty acids & Low-BMI"="#0F9D58")

env1_env2_prev_plt <- ggplot(env1_env2_prs_prev_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of PGS-Major depressive disorder") + ylab("Prevalence of disease (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) + ylim(0,20) + 
  scale_fill_manual(name="Group",values=model_colors)+
  scale_color_manual(name="Group",values=model_colors)
env1_env2_prev_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Figure6/O3FA_BMI_MajorDepression_double_interaction_PGSPrev_plt.pdf",height = 8,width = 8,dpi = 300,limitsize = FALSE)



final_plt <- ggarrange(low_high_o3fa_prev_plt, low_high_bmi_prev_plt, env1_env2_prev_plt,ncol = 3)
final_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Figure6/O3FA_BMI_MajorDepression_PGSPrev_plt.pdf",height = 8,width = 24,dpi = 300,limitsize = FALSE)
















