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


disease_status_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmerged_extractPhecodes.lst"))
disease_status_df$ID <- paste0(disease_status_df$EID,"_",disease_status_df$EID,sep="")
disease_status_df <- disease_status_df[c('ID','Status_574.1')]

prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/V2ALLmergedPRS_UKBB.tsv"))
prs_df <- prs_df[c('ID','PRS_finngen_R5_K11_CHOLELITH_PRS')]
colnames(prs_df) <- c('ID','PRS_ieu_a_1188')

covar_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/12_12_2022_UKBBCovariate_table.tsv"))
covar_df$ID <- paste0(covar_df$eid,"_",covar_df$eid,sep="")
covar_df$sex_corrected <- ifelse(covar_df$sex=="Male",1,0)
covar_df <- covar_df[c('ID','sex','age','age_square','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','sex_corrected')]

metabolite_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/Zscore_format_normalized_metabolite_table2_participant.tsv"))
metabolite_df$ID <- paste0(metabolite_df$eid,"_",metabolite_df$eid,sep="")
metabolite_df <- metabolite_df[c('ID','p23445_zscore')]
metabolite_df <- metabolite_df[!is.na(metabolite_df$p23445_zscore), ]

merge_df <- merge(covar_df,prs_df,by="ID")
merge_df <- merge(merge_df,disease_status_df,by="ID")
merge_df <- merge_df[merge_df$Status_574.1 != "Excluded", ]
merge_df$disease_status <- ifelse(merge_df$Status_574.1 == "Case", 1, 0)
MASTER_metabolite_df <- merge(merge_df,metabolite_df,by="ID")


# ------------------- Assess high vs. low DHA for Depression -----------------------
MASTER_metabolite_df$p23450_zscore_group <- ifelse(MASTER_metabolite_df$p23445_zscore >= mean(MASTER_metabolite_df$p23445_zscore),1,0) # 1 is HIGH DHA and 0 id LOW DHA

metabolite_gxe_glm = glm(disease_status~PRS_ieu_a_1188*p23450_zscore_group+PRS_ieu_a_1188+p23450_zscore_group+sex_corrected+age+age_square+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=MASTER_metabolite_df,family="binomial")
summary(metabolite_gxe_glm)

# Get low DHA
low_dha_df <- MASTER_metabolite_df[MASTER_metabolite_df$p23450_zscore_group == 0, ]
low_dha_df <- low_dha_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model1_prev_perc_df <- low_dha_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model1_prev_perc_df <- as.data.frame(model1_prev_perc_df)
model1_prev_perc_df$model <- "Low-Omega-6 fatty acids"
colnames(model1_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model1_prev_perc_df <- model1_prev_perc_df[c('PRS_percentile','Prev','model')]

# Get high DHA
high_dha_df <- MASTER_metabolite_df[MASTER_metabolite_df$p23450_zscore_group == 1, ]
high_dha_df <- high_dha_df %>% mutate(Percentile_PRS_ieu_a_1188 = ntile(PRS_ieu_a_1188, 100))

model2_prev_perc_df <- high_dha_df %>%
  group_by(Percentile_PRS_ieu_a_1188) %>%
  dplyr::summarize(Cases = sum(disease_status == 1), Controls = sum(disease_status == 0),Prev=(sum(disease_status == 1)/(sum(disease_status == 1)+sum(disease_status == 0))))

model2_prev_perc_df <- as.data.frame(model2_prev_perc_df)
model2_prev_perc_df$model <- "High-Omega-6 fatty acids"
colnames(model2_prev_perc_df) <- c('PRS_percentile','Cases','Controls','Prev','model')
model2_prev_perc_df <- model2_prev_perc_df[c('PRS_percentile','Prev','model')]

low_high_DHA_prs_prev_df <- rbind(model1_prev_perc_df,model2_prev_perc_df)

model_colors <- c("Low-Omega-6 fatty acids" = "#097969","High-Omega-6 fatty acids" = "#AFE1AF")

low_high_dha_prev_plt <- ggplot(low_high_DHA_prs_prev_df, aes(x=PRS_percentile, y=Prev*100,col = model,fill=model, group = model)) + geom_point(alpha=0.3) + 
  theme_bw() + theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  xlab("Percentile of PGS-Cholelithiasis") + ylab("Prevalence of disease (%)") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1) +
  scale_fill_manual(name="Group",values=model_colors)+ 
  scale_color_manual(name="Group",values=model_colors)
low_high_dha_prev_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/PRSxMetabolite/O6FA_Cholelithiasis_PGSPrev_plt.pdf",height = 8,width = 8,dpi = 300,limitsize = FALSE)
