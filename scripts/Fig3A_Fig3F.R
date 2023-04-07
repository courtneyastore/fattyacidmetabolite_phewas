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

map_phecode_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/phecode_definitions1.2.csv"))
map_phecode_df <- map_phecode_df[c("phecode","phenotype","category")]
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
map_phecode_df$category <- capFirst(map_phecode_df$category)

#------------------ Non-genetic O3FA--------------------------
o3fa_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/O3FAzscoreBMI_Interaction_logisticRegression_MRsupportdiseases.tsv"))
o3fa_bmi_genetic_df <- o3fa_bmi_genetic_df[c('phecode_id','zscoreBMI_zscoreO3FA_pval')]
colnames(o3fa_bmi_genetic_df) <- c("phecode","p")
o3fa_bmi_genetic_df$environment <- "BMI"
o3fa_bmi_genetic_df <- merge(o3fa_bmi_genetic_df,map_phecode_df,by = "phecode")
nrow(o3fa_bmi_genetic_df[o3fa_bmi_genetic_df$p < 0.05, ])

o3fa_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/O3FAzscoreWHR_Interaction_logisticRegression_MRsupportdiseases.tsv"))
o3fa_whr_genetic_df <- o3fa_whr_genetic_df[c('phecode_id','zscoreWHR_zscoreO3FA_pval')]
colnames(o3fa_whr_genetic_df) <- c("phecode","p")
o3fa_whr_genetic_df$environment <- "WHR"
o3fa_whr_genetic_df <- merge(o3fa_whr_genetic_df,map_phecode_df,by = "phecode")
nrow(o3fa_whr_genetic_df[o3fa_whr_genetic_df$p < 0.05, ])

o3fa_final_df <- rbind(o3fa_whr_genetic_df,o3fa_bmi_genetic_df)
o3fa_final_df$p_label <- ifelse(o3fa_final_df$p < 0.05, "*","")
o3fa_final_df$p_label <- ifelse(o3fa_final_df$p_label == "*" & o3fa_final_df$p < 0.001, "**",o3fa_final_df$p_label)
o3fa_final_df$p_label <- ifelse(o3fa_final_df$p_label == "**" & o3fa_final_df$p < 0.0001, "***",o3fa_final_df$p_label)

o3fa_plt <- ggplot(o3fa_final_df, aes(x=phenotype, y=environment)) + 
  geom_tile(aes(fill= -log10(p))) + 
  geom_text(aes(label=p_label),size=10) +
  xlab("") + ylab("") + ggtitle("Z-score-Omega-3 fatty acids x Anthropometric measurements") + 
  theme_classic() +  
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1),legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_fill_gradient2(low="white", high="blue4", guide="colorbar",limits=c(0,22.5)) + facet_grid(. ~ category, scales = "free_x", space='free')
o3fa_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric/O3FAZscore_anthropometric.pdf",plot=o3fa_plt,height=9,width = 38,dpi=300)

#------------------ Genetic O3FA--------------------------
pgs_o3fa_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/PGSO3FAzscoreBMI_Interaction_logisticRegression_MRsupportdiseases.tsv"))
pgs_o3fa_bmi_genetic_df <- pgs_o3fa_bmi_genetic_df[c('phecode_id','zscoreBMI_scalePRS_pval')]
colnames(pgs_o3fa_bmi_genetic_df) <- c("phecode","p")
pgs_o3fa_bmi_genetic_df$environment <- "BMI"
pgs_o3fa_bmi_genetic_df <- merge(pgs_o3fa_bmi_genetic_df,map_phecode_df,by = "phecode")
nrow(pgs_o3fa_bmi_genetic_df[pgs_o3fa_bmi_genetic_df$p < 0.05, ])

pgs_o3fa_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/PGSO3FAzscoreWHR_Interaction_logisticRegression_MRsupportdiseases.tsv"))
pgs_o3fa_whr_genetic_df <- pgs_o3fa_whr_genetic_df[c('phecode_id','zscoreWHR_scalePRS_pval')]
colnames(pgs_o3fa_whr_genetic_df) <- c("phecode","p")
pgs_o3fa_whr_genetic_df$environment <- "WHR"
pgs_o3fa_whr_genetic_df <- merge(pgs_o3fa_whr_genetic_df,map_phecode_df,by = "phecode")
nrow(pgs_o3fa_whr_genetic_df[pgs_o3fa_whr_genetic_df$p < 0.05, ])

pgs_o3fa_final_df <- rbind(pgs_o3fa_whr_genetic_df,pgs_o3fa_bmi_genetic_df)
pgs_o3fa_final_df$p_label <- ifelse(pgs_o3fa_final_df$p < 0.05, "*","")
pgs_o3fa_final_df$p_label <- ifelse(pgs_o3fa_final_df$p_label == "*" & pgs_o3fa_final_df$p < 0.001, "**",pgs_o3fa_final_df$p_label)
pgs_o3fa_final_df$p_label <- ifelse(pgs_o3fa_final_df$p_label == "**" & pgs_o3fa_final_df$p < 0.0001, "***",pgs_o3fa_final_df$p_label)

pgs_o3fa_plt <- ggplot(pgs_o3fa_final_df, aes(x=phenotype, y=environment)) + 
  geom_tile(aes(fill= -log10(p))) + 
  geom_text(aes(label=p_label),size=10) +
  xlab("") + ylab("") + ggtitle("PGS-Omega-3 fatty acids x Anthropometric measurements") + 
  theme_classic() +  
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1),legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_fill_gradient2(low="white", high="blue4", guide="colorbar",limits=c(0,22.5)) + facet_grid(. ~ category, scales = "free_x", space='free')
pgs_o3fa_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric/O3FAPGS_anthropometric.pdf",plot=pgs_o3fa_plt,height=9,width = 38,dpi=300)


#------------------ Non-genetic O6FA--------------------------
o6fa_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/O6FAzscoreBMI_Interaction_logisticRegression_MRsupportdiseases.tsv"))
o6fa_bmi_genetic_df <- o6fa_bmi_genetic_df[c('phecode_id','zscoreBMI_zscoreO6FA_pval')]
colnames(o6fa_bmi_genetic_df) <- c("phecode","p")
o6fa_bmi_genetic_df$environment <- "BMI"
o6fa_bmi_genetic_df <- merge(o6fa_bmi_genetic_df,map_phecode_df,by = "phecode")
nrow(o6fa_bmi_genetic_df[o6fa_bmi_genetic_df$p < 0.05, ])

o6fa_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/O6FAzscoreWHR_Interaction_logisticRegression_MRsupportdiseases.tsv"))
o6fa_whr_genetic_df <- o6fa_whr_genetic_df[c('phecode_id','zscoreWHR_zscoreO6FA_pval')]
colnames(o6fa_whr_genetic_df) <- c("phecode","p")
o6fa_whr_genetic_df$environment <- "WHR"
o6fa_whr_genetic_df <- merge(o6fa_whr_genetic_df,map_phecode_df,by = "phecode")
nrow(o6fa_whr_genetic_df[o6fa_whr_genetic_df$p < 0.05, ])


o6fa_final_df <- rbind(o6fa_whr_genetic_df,o6fa_bmi_genetic_df)
o6fa_final_df$p_label <- ifelse(o6fa_final_df$p < 0.05, "*","")
o6fa_final_df$p_label <- ifelse(o6fa_final_df$p_label == "*" & o6fa_final_df$p < 0.001, "**",o6fa_final_df$p_label)
o6fa_final_df$p_label <- ifelse(o6fa_final_df$p_label == "**" & o6fa_final_df$p < 0.0001, "***",o6fa_final_df$p_label)

o6fa_plt <- ggplot(o6fa_final_df, aes(x=phenotype, y=environment)) + 
  geom_tile(aes(fill= -log10(p))) + 
  geom_text(aes(label=p_label),size=10) +
  xlab("") + ylab("") + ggtitle("Z-score-Omega-6 fatty acids x Anthropometric measurements") + 
  theme_classic() +  
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1),legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_fill_gradient2(low="white", high="#097969", guide="colorbar",limits=c(0,22.5)) + facet_grid(. ~ category, scales = "free_x", space='free')
o6fa_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric/O6FAZscore_anthropometric.pdf",plot=o6fa_plt,height=9,width = 38,dpi=300)


#------------------ Genetic O6FA--------------------------
pgs_o6fa_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/PGSO6FAzscoreBMI_Interaction_logisticRegression_MRsupportdiseases.tsv"))
pgs_o6fa_bmi_genetic_df <- pgs_o6fa_bmi_genetic_df[c('phecode_id','zscoreBMI_scalePRS_pval')]
colnames(pgs_o6fa_bmi_genetic_df) <- c("phecode","p")
pgs_o6fa_bmi_genetic_df$environment <- "BMI"
pgs_o6fa_bmi_genetic_df <- merge(pgs_o6fa_bmi_genetic_df,map_phecode_df,by = "phecode")
nrow(pgs_o6fa_bmi_genetic_df[pgs_o6fa_bmi_genetic_df$p < 0.05, ])

pgs_o6fa_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/PGSO6FAzscoreWHR_Interaction_logisticRegression_MRsupportdiseases.tsv"))
pgs_o6fa_whr_genetic_df <- pgs_o6fa_whr_genetic_df[c('phecode_id','zscoreWHR_scalePRS_pval')]
colnames(pgs_o6fa_whr_genetic_df) <- c("phecode","p")
pgs_o6fa_whr_genetic_df$environment <- "WHR"
pgs_o6fa_whr_genetic_df <- merge(pgs_o6fa_whr_genetic_df,map_phecode_df,by = "phecode")
nrow(pgs_o6fa_whr_genetic_df[pgs_o6fa_whr_genetic_df$p < 0.05, ])


pgs_o6fa_final_df <- rbind(pgs_o6fa_whr_genetic_df,pgs_o6fa_bmi_genetic_df)
pgs_o6fa_final_df$p_label <- ifelse(pgs_o6fa_final_df$p < 0.05, "*","")
pgs_o6fa_final_df$p_label <- ifelse(pgs_o6fa_final_df$p_label == "*" & pgs_o6fa_final_df$p < 0.001, "**",pgs_o6fa_final_df$p_label)
pgs_o6fa_final_df$p_label <- ifelse(pgs_o6fa_final_df$p_label == "**" & pgs_o6fa_final_df$p < 0.0001, "***",pgs_o6fa_final_df$p_label)

pgs_o6fa_plt <- ggplot(pgs_o6fa_final_df, aes(x=phenotype, y=environment)) + 
  geom_tile(aes(fill= -log10(p))) + 
  geom_text(aes(label=p_label),size=10) +
  xlab("") + ylab("") + ggtitle("PGS-Omega-6 fatty acids x Anthropometric measurements") + 
  theme_classic() +  
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1),legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_fill_gradient2(low="white", high="#097969", guide="colorbar",limits=c(0,22.5)) + facet_grid(. ~ category, scales = "free_x", space='free')
pgs_o6fa_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric/O6FAPGS_anthropometric.pdf",plot=pgs_o6fa_plt,height=9,width = 38,dpi=300)


#------------------ Non-genetic DHA--------------------------
dha_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/DHAzscoreBMI_Interaction_logisticRegression_MRsupportdiseases.tsv"))
dha_bmi_genetic_df <- dha_bmi_genetic_df[c('phecode_id','zscoreBMI_zscoreDHA_pval')]
colnames(dha_bmi_genetic_df) <- c("phecode","p")
dha_bmi_genetic_df$environment <- "BMI"
dha_bmi_genetic_df <- merge(dha_bmi_genetic_df,map_phecode_df,by = "phecode")
nrow(dha_bmi_genetic_df[dha_bmi_genetic_df$p < 0.05, ])

dha_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/DHAzscoreWHR_Interaction_logisticRegression_MRsupportdiseases.tsv"))
dha_whr_genetic_df <- dha_whr_genetic_df[c('phecode_id','zscoreWHR_zscoreDHA_pval')]
colnames(dha_whr_genetic_df) <- c("phecode","p")
dha_whr_genetic_df$environment <- "WHR"
dha_whr_genetic_df <- merge(dha_whr_genetic_df,map_phecode_df,by = "phecode")
nrow(dha_whr_genetic_df[dha_whr_genetic_df$p < 0.05, ])

dha_final_df <- rbind(dha_whr_genetic_df,dha_bmi_genetic_df)
dha_final_df$p_label <- ifelse(dha_final_df$p < 0.05, "*","")
dha_final_df$p_label <- ifelse(dha_final_df$p_label == "*" & dha_final_df$p < 0.001, "**",dha_final_df$p_label)
dha_final_df$p_label <- ifelse(dha_final_df$p_label == "**" & dha_final_df$p < 0.0001, "***",dha_final_df$p_label)

dha_plt <- ggplot(dha_final_df, aes(x=phenotype, y=environment)) + 
  geom_tile(aes(fill= -log10(p))) + 
  geom_text(aes(label=p_label),size=10) +
  xlab("") + ylab("") + ggtitle("Z-score-Docosahexaenoic acid x Anthropometric measurements") + 
  theme_classic() +  
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1),legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_fill_gradient2(low="white", high="#ff781f", guide="colorbar",limits=c(0,22.5)) + facet_grid(. ~ category, scales = "free_x", space='free')
dha_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric/DHAZscore_anthropometric.pdf",plot=dha_plt,height=9,width = 65,dpi=300,limitsize = FALSE)

#------------------ Genetic DHA--------------------------
pgs_dha_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/PGSDHAzscoreBMI_Interaction_logisticRegression_MRsupportdiseases.tsv"))
pgs_dha_bmi_genetic_df <- pgs_dha_bmi_genetic_df[c('phecode_id','zscoreBMI_scalePRS_pval')]
colnames(pgs_dha_bmi_genetic_df) <- c("phecode","p")
pgs_dha_bmi_genetic_df$environment <- "BMI"
pgs_dha_bmi_genetic_df <- merge(pgs_dha_bmi_genetic_df,map_phecode_df,by = "phecode")
nrow(pgs_dha_bmi_genetic_df[pgs_dha_bmi_genetic_df$p < 0.05, ])


pgs_dha_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/BMIWHR_results/PGSDHAzscoreWHR_Interaction_logisticRegression_MRsupportdiseases.tsv"))
pgs_dha_whr_genetic_df <- pgs_dha_whr_genetic_df[c('phecode_id','zscoreWHR_scalePRS_pval')]
colnames(pgs_dha_whr_genetic_df) <- c("phecode","p")
pgs_dha_whr_genetic_df$environment <- "WHR"
pgs_dha_whr_genetic_df <- merge(pgs_dha_whr_genetic_df,map_phecode_df,by = "phecode")
nrow(pgs_dha_whr_genetic_df[pgs_dha_whr_genetic_df$p < 0.05, ])


pgs_dha_final_df <- rbind(pgs_dha_whr_genetic_df,pgs_dha_bmi_genetic_df)
#pgs_dha_final_df$p_label <- ifelse(pgs_dha_final_df$p < 0.05, round(-log10(pgs_dha_final_df$p),2),"")
pgs_dha_final_df$p_label <- ifelse(pgs_dha_final_df$p < 0.05, "*","")
pgs_dha_final_df$p_label <- ifelse(pgs_dha_final_df$p_label == "*" & pgs_dha_final_df$p < 0.001, "**",pgs_dha_final_df$p_label)
pgs_dha_final_df$p_label <- ifelse(pgs_dha_final_df$p_label == "**" & pgs_dha_final_df$p < 0.0001, "***",pgs_dha_final_df$p_label)


pgs_dha_plt <- ggplot(pgs_dha_final_df, aes(x=phenotype, y=environment)) + 
  geom_tile(aes(fill= -log10(p))) + 
  geom_text(aes(label=p_label),size=10) +
  xlab("") + ylab("") + ggtitle("PGS-Docosahexaenoic acid x Anthropometric measurements") + 
  theme_classic() +  
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1),legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_fill_gradient2(low="white", high="#ff781f", guide="colorbar",limits=c(0,22.5)) + facet_grid(. ~ category, scales = "free_x", space='free')
pgs_dha_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric/DHAPGS_anthropometric.pdf",plot=pgs_dha_plt,height=9,width = 65,dpi=300,limitsize = FALSE)














#o6fa_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmerged_extractPhecodes.lst"))
#dha_bmi_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmerged_extractPhecodes.lst"))

#o3fa_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmerged_extractPhecodes.lst"))
#o6fa_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmerged_extractPhecodes.lst"))
#dha_whr_genetic_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ALLmerged_extractPhecodes.lst"))

ggplot(o3fa_final_df, aes(y=phenotype, x=environment, fill= -log10(p))) + 
  xlab("") + ylab("") +  ggtitle("Omega-3 fatty acid levels x Anthropometric measurements") + 
  theme_classic() + 
  theme(legend.position="bottom",text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  geom_tile() + scale_fill_gradient2(low="white", high="blue4", guide="colorbar")

