library(ggplot2)
library(forcats)
library(tidyverse)
library(dplyr)
library(data.table)
library("ggpubr")
library(plyr)
library(stringr)
library(rsq)

disease_status_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/V4_ALLmerged_extractPhecodes.lst"))
disease_status_df$ID <- paste0(disease_status_df$EID,"_",disease_status_df$EID,sep="")

met_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/Zscore_format_normalized_metabolite_table2_participant.tsv"))
met_df$ID <- paste0(met_df$eid,"_",met_df$eid,sep="")
met_df <- met_df[c('ID','p23450_zscore')]

covar_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/12_12_2022_UKBBCovariate_table.tsv"))
covar_df$ID <- paste0(covar_df$eid,"_",covar_df$eid,sep="")
covar_df$sex_corrected <- ifelse(covar_df$sex=="Male",1,0)

merge_df <- merge(covar_df,met_df,by="ID")
merge_df <- merge(merge_df,disease_status_df,by="ID")
merge_df <- merge_df[!is.na(merge_df$p23450_zscore), ]

merge_df$disease_status <- ifelse(merge_df$Status_296.22 == "Case",1,0)

# Run logistic regression model
glm(disease_status ~ p23450_zscore + age + age_square + sex_corrected + PC1 + PC2 + PC3 + PC4 + PC5, data = merge_df)

# Calculate the percentiles for the PRS column
merge_df <- merge_df %>% mutate(Percentile_PRS_met_d_dha= ntile(p23450_zscore, 100))

model1_prev_perc_df <- merge_df %>%
  group_by(Percentile_PRS_met_d_dha) %>%
  dplyr::summarize(Cases = sum(Status_296.22 == "Case"), Controls = sum(Status_296.22 == "Control"),Prev=(sum(Status_296.22 == "Case")/(sum(Status_296.22 == "Case")+sum(Status_296.22 == "Control"))))

model1_prev_perc_df <- as.data.frame(model1_prev_perc_df)
model1_prev_perc_df$model <- "Major depressive disorder"
colnames(model1_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')

# Generate figure
prev_plt <- ggplot(model1_prev_perc_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + theme(text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of Docosahexaenoic acid Z-score") + ylab("Prevalence of disease (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) + ylim(0,11.5)
prev_plt
ggsave("DHAZscorePerc_PrevDepression_plt.pdf",height = 6,width = 12,dpi = 300,limitsize = FALSE)










