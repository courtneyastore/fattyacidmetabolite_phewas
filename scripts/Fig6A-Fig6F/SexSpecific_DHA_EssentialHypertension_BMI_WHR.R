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
prs_df <- prs_df[c('ID','PRS_met-d-dha_PRS')]
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

# Obese female BMI
obese_bmi_female_df <- merge_df[(merge_df$bmi_group == 1) & (merge_df$sex_corrected == 0), ]
obese_bmi_female_df <- obese_bmi_female_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model1_prev_perc_df <- obese_bmi_female_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model1_prev_perc_df <- as.data.frame(model1_prev_perc_df)
model1_prev_perc_df$model <- "Sex-Female & BMI-Obese"
colnames(model1_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model1_prev_perc_df <- model1_prev_perc_df[c('PRS_percentile','Prev','model')]

# Obese male BMI
obese_bmi_male_df <- merge_df[(merge_df$bmi_group == 1) & (merge_df$sex_corrected == 1), ]
obese_bmi_male_df <- obese_bmi_male_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model2_prev_perc_df <- obese_bmi_male_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model2_prev_perc_df <- as.data.frame(model2_prev_perc_df)
model2_prev_perc_df$model <- "Sex-Male & BMI-Obese"
colnames(model2_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model2_prev_perc_df <- model2_prev_perc_df[c('PRS_percentile','Prev','model')]

# Non-obese female BMI
non_obese_bmi_female_df <- merge_df[(merge_df$bmi_group == 0) & (merge_df$sex_corrected == 0), ]
non_obese_bmi_female_df <- non_obese_bmi_female_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model3_prev_perc_df <- non_obese_bmi_female_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model3_prev_perc_df <- as.data.frame(model3_prev_perc_df)
model3_prev_perc_df$model <- "Sex-Female & BMI-Non-obese"
colnames(model3_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model3_prev_perc_df <- model3_prev_perc_df[c('PRS_percentile','Prev','model')]

# Non-obese female BMI
non_obese_bmi_male_df <- merge_df[(merge_df$bmi_group == 0) & (merge_df$sex_corrected == 1), ]
non_obese_bmi_male_df <- non_obese_bmi_male_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model4_prev_perc_df <- non_obese_bmi_male_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model4_prev_perc_df <- as.data.frame(model4_prev_perc_df)
model4_prev_perc_df$model <- "Sex-Male & BMI-Non-obese"
colnames(model4_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model4_prev_perc_df <- model4_prev_perc_df[c('PRS_percentile','Prev','model')]


# Obese female WHR
obese_whr_female_df <- merge_df[(merge_df$whr_group == 1) & (merge_df$sex_corrected == 0), ]
obese_whr_female_df <- obese_whr_female_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model5_prev_perc_df <- obese_whr_female_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model5_prev_perc_df <- as.data.frame(model5_prev_perc_df)
model5_prev_perc_df$model <- "Sex-Female & WHR-Obese"
colnames(model5_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model5_prev_perc_df <- model5_prev_perc_df[c('PRS_percentile','Prev','model')]

# Obese male WHR
obese_whr_male_df <- merge_df[(merge_df$whr_group == 1) & (merge_df$sex_corrected == 1), ]
obese_whr_male_df <- obese_whr_male_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model6_prev_perc_df <- obese_whr_male_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model6_prev_perc_df <- as.data.frame(model6_prev_perc_df)
model6_prev_perc_df$model <- "Sex-Male & WHR-Obese"
colnames(model6_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model6_prev_perc_df <- model6_prev_perc_df[c('PRS_percentile','Prev','model')]

# Non-obese female WHR
non_obese_whr_female_df <- merge_df[(merge_df$whr_group == 0) & (merge_df$sex_corrected == 0), ]
non_obese_whr_female_df <- non_obese_whr_female_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model7_prev_perc_df <- non_obese_whr_female_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model7_prev_perc_df <- as.data.frame(model7_prev_perc_df)
model7_prev_perc_df$model <- "Sex-Female & WHR-Non-obese"
colnames(model7_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model7_prev_perc_df <- model7_prev_perc_df[c('PRS_percentile','Prev','model')]

# Non-obese male WHR
non_obese_whr_male_df <- merge_df[(merge_df$whr_group == 0) & (merge_df$sex_corrected == 1), ]
non_obese_whr_male_df <- non_obese_whr_male_df %>% mutate(Percentile_PRS = ntile(PRS, 100))
model8_prev_perc_df <- non_obese_whr_male_df %>%
  group_by(Percentile_PRS) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))
model8_prev_perc_df <- as.data.frame(model8_prev_perc_df)
model8_prev_perc_df$model <- "Sex-Male & WHR-Non-obese"
colnames(model8_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model8_prev_perc_df <- model8_prev_perc_df[c('PRS_percentile','Prev','model')]

final_df <- rbind(model1_prev_perc_df,model2_prev_perc_df,model3_prev_perc_df,model4_prev_perc_df)

bmi_plt <- ggplot(final_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + ggtitle("BMI") + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of PGS-Docosahexaenoic acid") + ylab("Prevalence of Essential hypertension (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1)
bmi_plt

final_df <- rbind(model5_prev_perc_df,model6_prev_perc_df,model7_prev_perc_df,model8_prev_perc_df)
whr_plt <- ggplot(final_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + ggtitle("WHR") + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of PGS-Docosahexaenoic acid") + ylab("Prevalence of Essential hypertension (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1)
whr_plt

# basic models
test_whr <- glm(data=merge_df,disease_status ~ WaistCircumference_median_to_HipCircumference_median_ratio + PRS + WaistCircumference_median_to_HipCircumference_median_ratio*PRS + age + age_square + sex_corrected + PC1 + PC2 + PC3 + PC4 + PC5, family = "binomial")
summary(test_whr)
test_bmi <- glm(data=merge_df,disease_status ~ BMI_median + PRS + BMI_median*PRS + age + age_square + sex_corrected + PC1 + PC2 + PC3 + PC4 + PC5, family = "binomial")
summary(test_bmi)

# no sex
test_nosex_whr <- glm(data=merge_df,disease_status ~ WaistCircumference_median_to_HipCircumference_median_ratio + PRS + WaistCircumference_median_to_HipCircumference_median_ratio*PRS + age + age_square + PC1 + PC2 + PC3 + PC4 + PC5, family = "binomial")
summary(test_nosex_whr)
test_nosex_bmi <- glm(data=merge_df,disease_status ~ BMI_median + PRS + BMI_median*PRS + age + age_square + PC1 + PC2 + PC3 + PC4 + PC5, family = "binomial")
summary(test_nosex_bmi)

# no covariates 
test_whr_nocovar <- glm(data=merge_df,disease_status ~ WaistCircumference_median_to_HipCircumference_median_ratio + PRS + WaistCircumference_median_to_HipCircumference_median_ratio*PRS, family = "binomial")
summary(test_whr_nocovar)
test_bmi_nocovar <- glm(data=merge_df,disease_status ~ BMI_median + PRS + BMI_median*PRS, family = "binomial")
summary(test_bmi_nocovar)


# Obese/non-obese status
test_whr_group <- glm(data=merge_df,disease_status ~ whr_group + PRS + whr_group*PRS + age + age_square + sex_corrected + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial")
summary(test_whr_group)

test_bmi_group <- glm(data=merge_df,disease_status ~ bmi_group + PRS + bmi_group*PRS + age + age_square + sex_corrected + PC1 + PC2 + PC3 + PC4 + PC5+ PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial")
summary(test_bmi_group)
