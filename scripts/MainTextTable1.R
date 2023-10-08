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

o3fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST8_suggestive_O3FA.tsv"))

o6fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST9_suggestive_O6FA.tsv"))

dha_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST10_suggestive_DHA.tsv"))

mr_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/MR_Batch_Analysis/results/MR_pairs_fatty_acid_SuggestedAssociations.tsv"))

# MR
mr_o3fa_df <- mr_df[(mr_df$id.exposure == "met-d-Omega_3") & (mr_df$method == "Inverse variance weighted"), ]
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
colnames(mr_o3fa_df) <- c('metabolite','disease','disease_group','MR_lo_ci','MR_hi_ci','MR_or','MR_pval')
mr_o3fa_df <- mr_o3fa_df[(mr_o3fa_df$MR_pval < 0.05) & (mr_o3fa_df$MR_or < 1), ]
mr_o3fa_df <- mr_o3fa_df[c('disease','MR_lo_ci','MR_hi_ci','MR_or','MR_pval')]
mr_o3fa_df <- merge(mr_o3fa_df,o3fa_df,by="disease")

# Number of O3FA MR disease associations
nrow(mr_o3fa_df)

# MR O6fa
mr_o6fa_df <- mr_df[(mr_df$id.exposure == "met-d-Omega_6") & (mr_df$method == "Inverse variance weighted"), ]
mr_o6fa_df <- mr_o6fa_df[c('exposure','id.outcome','or_lci95','or_uci95','or','pval')]
mr_o6fa_df$metabolte <- gsub("\\|.*", "", mr_o6fa_df$exposure)
mr_o6fa_df$metabolte <- str_trim(mr_o6fa_df$metabolte)
mr_o6fa_df <- mr_o6fa_df[c('metabolte','id.outcome','or_lci95','or_uci95','or','pval')]
colnames(mr_o6fa_df) <- c('metabolite','GWAS_id','lo_ci','hi_ci','or','p')
mr_o6fa_df <- merge(map_gwasid_df,mr_o6fa_df,by="GWAS_id")
mr_o6fa_df <- unique(mr_o6fa_df)
colnames(mr_o6fa_df) <- c('GWAS_id','phenotype','metabolite','lo_ci','hi_ci','or','p')
mr_o6fa_df <- merge(mr_o6fa_df,phecode_map_df,by="phenotype")
mr_o6fa_df <- mr_o6fa_df[c('metabolite','phenotype','category','lo_ci','hi_ci','or','p')]
colnames(mr_o6fa_df) <- c('metabolite','disease','disease_group','MR_lo_ci','MR_hi_ci','MR_or','MR_pval')
mr_o6fa_df <- mr_o6fa_df[(mr_o6fa_df$MR_pval < 0.05) & (mr_o6fa_df$MR_or < 1), ]
mr_o6fa_df <- mr_o6fa_df[c('disease','MR_lo_ci','MR_hi_ci','MR_or','MR_pval')]
mr_o6fa_df <- merge(mr_o6fa_df,o6fa_df,by="disease")

# Number of O3FA MR disease associations
nrow(mr_o6fa_df)

# MR DHA
mr_dha_df <- mr_df[(mr_df$id.exposure == "met-d-DHA") & (mr_df$method == "Inverse variance weighted"), ]
mr_dha_df <- mr_dha_df[c('exposure','id.outcome','or_lci95','or_uci95','or','pval')]
mr_dha_df$metabolte <- gsub("\\|.*", "", mr_dha_df$exposure)
mr_dha_df$metabolte <- str_trim(mr_dha_df$metabolte)
mr_dha_df <- mr_dha_df[c('metabolte','id.outcome','or_lci95','or_uci95','or','pval')]
colnames(mr_dha_df) <- c('metabolite','GWAS_id','lo_ci','hi_ci','or','p')
mr_dha_df <- merge(map_gwasid_df,mr_dha_df,by="GWAS_id")
mr_dha_df <- unique(mr_dha_df)
colnames(mr_dha_df) <- c('GWAS_id','phenotype','metabolite','lo_ci','hi_ci','or','p')
mr_dha_df <- merge(mr_dha_df,phecode_map_df,by="phenotype")
mr_dha_df <- mr_dha_df[c('metabolite','phenotype','category','lo_ci','hi_ci','or','p')]
colnames(mr_dha_df) <- c('metabolite','disease','disease_group','MR_lo_ci','MR_hi_ci','MR_or','MR_pval')
mr_dha_df <- mr_dha_df[(mr_dha_df$MR_pval < 0.05) & (mr_dha_df$MR_or < 1), ]
mr_dha_df <- mr_dha_df[c('disease','MR_lo_ci','MR_hi_ci','MR_or','MR_pval')]
mr_dha_df <- merge(mr_dha_df,dha_df,by="disease")

# Number of O3FA MR disease associations
nrow(mr_dha_df)

table1_df <- rbind(mr_o3fa_df,mr_o6fa_df,mr_dha_df)
write.table(table1_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/MainText_table1.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


