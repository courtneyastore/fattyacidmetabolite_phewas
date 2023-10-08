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
library(viridis)
library(hrbrthemes)
library(circlize)
circos.par("track.height" = 0.4)
options(ggrepel.max.overlaps = Inf)

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

p_thresh = 0.05/(177*3)

map_gwasid_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/All_SuggestiveDiseases.txt"))
phecode_map_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/phecode_definitions1.2.csv"))
phecode_map_df$phecode <- as.character(phecode_map_df$phecode)
phecode_map_df$category <- capFirst(phecode_map_df$category)


mr_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/MR_Batch_Analysis/results/MR_pairs_fatty_acid_SuggestedAssociations.tsv"))
map_gwasid_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/All_SuggestiveDiseases.txt"))

# MR
mr_o3fa_df <- mr_df[mr_df$method == "Inverse variance weighted", ]
mr_o3fa_df <- mr_o3fa_df[c('exposure','id.outcome','or_lci95','or_uci95','or','pval')]
mr_o3fa_df$metabolte <- gsub("\\|.*", "", mr_o3fa_df$exposure)
mr_o3fa_df$metabolte <- str_trim(mr_o3fa_df$metabolte)
mr_o3fa_df <- mr_o3fa_df[c('metabolte','id.outcome','or_lci95','or_uci95','or','pval')]
colnames(mr_o3fa_df) <- c('metabolite','GWAS_id','lo_ci','hi_ci','or','p')
mr_o3fa_df <- merge(map_gwasid_df,mr_o3fa_df,by="GWAS_id")
mr_o3fa_df <- unique(mr_o3fa_df)
colnames(mr_o3fa_df) <- c('GWAS_id','phenotype','metabolite','lo_ci','hi_ci','or','p')
mr_o3fa_df <- merge(mr_o3fa_df,phecode_map_df,by="phenotype")
mr_o3fa_df <- mr_o3fa_df[c('metabolite','phenotype','category','lo_ci','hi_ci','or','p')]
colnames(mr_o3fa_df) <- c('metabolite','disease','category','lo_ci','hi_ci','or','pval')

mr_o3fa_df <- mr_o3fa_df[(mr_o3fa_df$p < 0.05), ]
mr_o3fa_df$metabolite <- ifelse(mr_o3fa_df$metabolite == "Omega-3 fatty acids","O3FA",mr_o3fa_df$metabolite)
mr_o3fa_df$metabolite <- ifelse(mr_o3fa_df$metabolite == "Omega-6 fatty acids","O6FA",mr_o3fa_df$metabolite)
mr_o3fa_df$metabolite <- ifelse(mr_o3fa_df$metabolite == "Docosahexaenoic acid","DHA",mr_o3fa_df$metabolite)

mr_o3fa_df$lo_ci <- format(mr_o3fa_df$lo_ci, scientific = TRUE, digits = 2)
mr_o3fa_df$hi_ci <- format(mr_o3fa_df$hi_ci, scientific = TRUE, digits = 2)
mr_o3fa_df$or <- format(mr_o3fa_df$or, scientific = TRUE, digits = 2)
mr_o3fa_df$pval <- format(mr_o3fa_df$pval, scientific = TRUE, digits = 2)

write.table(mr_o3fa_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/PUFA_canalization_RshinyApp/PUFA_MR_table.tsv",quote=TRUE,sep="\t",row.names=FALSE,col.names = TRUE)
