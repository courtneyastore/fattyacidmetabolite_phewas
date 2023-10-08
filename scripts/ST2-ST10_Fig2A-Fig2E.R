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

p_thresh = 0.05/(1326*6)

phecode_map_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/phecode_definitions1.2.csv"))
phecode_map_df$phecode <- as.character(phecode_map_df$phecode)
phecode_map_df$category <- capFirst(phecode_map_df$category)

# Read and format O3FA genetic results
o3fa_prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/result_files/MASTER_01_16_2023_Phecodes_PRS_O3FA_age_agesquare_sex_10PCS_glmresults.tsv"))
o3fa_prs_df$phecode <- as.character(o3fa_prs_df$phecode)
o3fa_prs_df$label <- "Omega-3 fatty acids"
o3fa_prs_df <- merge(o3fa_prs_df,phecode_map_df,by="phecode")
o3fa_prs_df <- o3fa_prs_df[o3fa_prs_df$category != "NULL", ]
o3fa_prs_df <- o3fa_prs_df[!(is.na(o3fa_prs_df$category)), ]
o3fa_prs_df <- o3fa_prs_df[o3fa_prs_df$category != "", ]
o3fa_prs_df$significance <- ifelse(o3fa_prs_df$pval< p_thresh, o3fa_prs_df$phenotype, "")
o3fa_prs_df$risk_prot <- ifelse(o3fa_prs_df$or > 1, "Risk", "Protective")

# Read and format O6FA genetic results
o6fa_prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/result_files/MASTER_01_16_2023_Phecodes_PRS_O6FA_age_agesquare_sex_10PCS_glmresults.tsv"))
o6fa_prs_df$phecode <- as.character(o6fa_prs_df$phecode)
o6fa_prs_df$label <- "Omega-6 fatty acids"
o6fa_prs_df <- merge(o6fa_prs_df,phecode_map_df,by="phecode")
o6fa_prs_df <- o6fa_prs_df[o6fa_prs_df$category != "NULL", ]
o6fa_prs_df <- o6fa_prs_df[!(is.na(o6fa_prs_df$category)), ]
o6fa_prs_df <- o6fa_prs_df[o6fa_prs_df$category != "", ]
o6fa_prs_df$significance <- ifelse(o6fa_prs_df$pval< p_thresh, o6fa_prs_df$phenotype, "")
o6fa_prs_df$risk_prot <- ifelse(o6fa_prs_df$or > 1, "Risk", "Protective")

# Read and format DHA genetic results
dha_prs_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/result_files/MASTER_01_16_2023_Phecodes_PRS_DHA_age_agesquare_sex_10PCS_glmresults.tsv"))
dha_prs_df$phecode <- as.character(dha_prs_df$phecode)
dha_prs_df$label <- "Docosahexaenoic acid"
dha_prs_df <- merge(dha_prs_df,phecode_map_df,by="phecode")
dha_prs_df <- dha_prs_df[dha_prs_df$category != "NULL", ]
dha_prs_df <- dha_prs_df[!(is.na(dha_prs_df$category)), ]
dha_prs_df <- dha_prs_df[dha_prs_df$category != "", ]
dha_prs_df$significance <- ifelse(dha_prs_df$pval< p_thresh, dha_prs_df$phenotype, "")
dha_prs_df$risk_prot <- ifelse(dha_prs_df$or > 1, "Risk", "Protective")

# Read and format O3FA results
o3fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/result_files/MASTER_01_16_2023_Phecodes_Zscore_O3FA_age_agesquare_sex_10PCS_glmresults.tsv"))
o3fa_df$phecode <- as.character(o3fa_df$phecode)
o3fa_df$label <- "Omega-3 fatty acids"
o3fa_df <- merge(o3fa_df,phecode_map_df,by="phecode")
o3fa_df <- o3fa_df[o3fa_df$category != "NULL", ]
o3fa_df <- o3fa_df[!(is.na(o3fa_df$category)), ]
o3fa_df <- o3fa_df[o3fa_df$category != "", ]
o3fa_df$significance <- ifelse(o3fa_df$pval< p_thresh, o3fa_df$phenotype, "")
o3fa_df$risk_prot <- ifelse(o3fa_df$or > 1, "Risk", "Protective")

# Read and format O6FA results
o6fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/result_files/MASTER_01_16_2023_Phecodes_Zscore_O6FA_age_agesquare_sex_10PCS_glmresults.tsv"))
o6fa_df$phecode <- as.character(o6fa_df$phecode)
o6fa_df$label <- "Omega-6 fatty acids"
o6fa_df <- merge(o6fa_df,phecode_map_df,by="phecode")
o6fa_df <- o6fa_df[o6fa_df$category != "NULL", ]
o6fa_df <- o6fa_df[!(is.na(o6fa_df$category)), ]
o6fa_df <- o6fa_df[o6fa_df$category != "", ]
o6fa_df$significance <- ifelse(o6fa_df$pval< p_thresh, o3fa_df$phenotype, "")
o6fa_df$risk_prot <- ifelse(o6fa_df$or > 1, "Risk", "Protective")

# Read and format DHA results
dha_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/result_files/MASTER_01_16_2023_Phecodes_Zscore_DHA_age_agesquare_sex_10PCS_glmresults.tsv"))
dha_df$phecode <- as.character(dha_df$phecode)
dha_df$label <- "Docosahexaenoic acid"
dha_df <- merge(dha_df,phecode_map_df,by="phecode")
dha_df <- dha_df[dha_df$category != "NULL", ]
dha_df <- dha_df[!(is.na(dha_df$category)), ]
dha_df <- dha_df[dha_df$category != "", ]
dha_df$significance <- ifelse(dha_df$pval< p_thresh, dha_df$phenotype, "")
dha_df$risk_prot <- ifelse(dha_df$or > 1, "Risk", "Protective")

#-------------------------------------------------------------------------

#--------------------Non-genetic association numbers----------------------
#--------------------------------O3FA-------------------------------------
# Number of total non-genetic associations O3FA
nrow(o3fa_df[o3fa_df$significance != "", ])

# Number of protective non-genetic associations O3FA
nrow(o3fa_df[(o3fa_df$significance != "") & (o3fa_df$risk_prot == "Protective"), ])

# Unique disease groups for protective associations O3FA
unique(o3fa_df[(o3fa_df$significance != "") & (o3fa_df$risk_prot == "Protective"), ]$category)

# Number of protective non-genetic associations O3FA
nrow(o3fa_df[(o3fa_df$significance != "") & (o3fa_df$risk_prot == "Risk"), ])

# Unique disease groups for risk associations O3FA
unique(o3fa_df[(o3fa_df$significance != "") & (o3fa_df$risk_prot == "Risk"), ]$category)

# Get significant O6FA associations
significant_o3fa_df <- o3fa_df[(o3fa_df$significance != "") & (o3fa_df$risk_prot == "Protective"), ]
significant_o3fa_df$label <- "Omega-3 fatty acids"
st1_df <- significant_o3fa_df[c("label","phenotype","category","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")]
colnames(st1_df) <- c("metabolite","disease","disease_group","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")
write.table(st1_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/V2_ST2_nongenetic_O3FA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


# Make table of DHA significant association category frequencies
category_significant_o3fa_df <- count(significant_o3fa_df,'category')
category_significant_o3fa_df$from <- "Omega-3 fatty acids"

#--------------------------------O6FA-------------------------------------
# Number of total non-genetic associations O6FA
nrow(o6fa_df[o6fa_df$significance != "", ])

# Number of protective non-genetic associations O6FA
nrow(o6fa_df[(o6fa_df$significance != "") & (o6fa_df$risk_prot == "Protective"), ])

# Unique disease groups for protective associations O6FA
unique(o6fa_df[(o6fa_df$significance != "") & (o6fa_df$risk_prot == "Protective"), ]$category)

# Number of protective non-genetic associations O6FA
nrow(o6fa_df[(o6fa_df$significance != "") & (o6fa_df$risk_prot == "Risk"), ])

# Unique disease groups for risk associations O6FA
unique(o6fa_df[(o6fa_df$significance != "") & (o6fa_df$risk_prot == "Risk"), ]$category)

# Get significant O6FA associations
significant_o6fa_df <- o6fa_df[(o6fa_df$significance != "") & (o6fa_df$risk_prot == "Protective"), ]
st2_df <- significant_o6fa_df[c("label","phenotype","category","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")]
colnames(st2_df) <- c("metabolite","disease","disease_group","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")
write.table(st2_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST2_nongenetic_O6FA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


# Make table of DHA significant association category frequencies
category_significant_o6fa_df <- count(significant_o6fa_df,'category')
category_significant_o6fa_df$from <- "Omega-6 fatty acids"


#--------------------------------DHA-------------------------------------
# Number of total non-genetic associations DHA
nrow(dha_df[dha_df$significance != "", ])

# Number of protective non-genetic associations DHA
nrow(dha_df[(dha_df$significance != "") & (dha_df$risk_prot == "Protective"), ])

# Unique disease groups for protective associations DHA
unique(dha_df[(dha_df$significance != "") & (dha_df$risk_prot == "Protective"), ]$category)

# Number of protective non-genetic associations DHA
nrow(dha_df[(dha_df$significance != "") & (dha_df$risk_prot == "Risk"), ])

# Unique disease groups for risk associations DHA
unique(dha_df[(dha_df$significance != "") & (dha_df$risk_prot == "Risk"), ]$category)

# Get significant DHA associations
significant_dha_df <- dha_df[(dha_df$significance != "") & (dha_df$risk_prot == "Protective"), ]
st3_df <- significant_dha_df[c("label","phenotype","category","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")]
colnames(st3_df) <- c("metabolite","disease","disease_group","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")
write.table(st3_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST3_nongenetic_DHA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


# Make table of DHA significant association category frequencies
category_significant_dha_df <- count(significant_dha_df,'category')
category_significant_dha_df$from <- "Docosahexaenoic acid"

# ------------------------Combine O3FA, O6FA, and DHA category counts for Circos plot----------------------------------
non_genetic_significant_protective_catcount_df <- rbind(category_significant_dha_df,category_significant_o6fa_df,category_significant_o3fa_df)
non_genetic_significant_protective_catcount_df <- non_genetic_significant_protective_catcount_df[c('from','category','freq')]
colnames(non_genetic_significant_protective_catcount_df) <- c("from","to","value")

grid.col = c('Omega-3 fatty acids'='blue4','Omega-6 fatty acids'='#097969','Docosahexaenoic acid'='#ff781f',
             'Digestive'="#5A5A5A",'Dermatologic'='Grey','Circulatory system'='#5A5A5A','Endocrine/metabolic'="Grey",
             'Genitourinary'="#5A5A5A","Hematopoietic"="Grey","Infectious diseases"="#5A5A5A","Injuries & poisonings"="Grey",
             "Mental disorders"="#5A5A5A","Musculoskeletal"="Grey","Neoplasms"="#5A5A5A","Neurological"="Grey",
             "Respiratory"="#5A5A5A","Sense organs"="Grey","Symptoms"="#5A5A5A")
order = c("Docosahexaenoic acid","Omega-6 fatty acids","Omega-3 fatty acids","Digestive", "Dermatologic", "Circulatory system", "Endocrine/metabolic", "Genitourinary", "Hematopoietic", 
          "Infectious diseases", "Injuries & poisonings", "Mental disorders","Musculoskeletal","Neoplasms",
          "Neurological","Respiratory","Sense organs","Symptoms")

pdf("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/NonGeneticAssociations_O3FA_O6FA_DHA_Protective_ChordPlt.pdf",height = 20,width = 20)
chordDiagram(non_genetic_significant_protective_catcount_df,grid.col = grid.col, order=order)
dev.off()

#-------------------Venn Diagram for non-genetic associations-----------------------------
non_genetic_sig_o3fa_diseases <- as.list(significant_o3fa_df$phenotype)
non_genetic_sig_o6fa_diseases <- as.list(significant_o6fa_df$phenotype)
non_genetic_sig_dha_diseases <- as.list(significant_dha_df$phenotype)

non_genetic_assoc <- list(
  "Omega-3 fatty acids" = non_genetic_sig_o3fa_diseases, 
  "Omega-6 fatty acids" = non_genetic_sig_o6fa_diseases, 
  "Docosahexaenoic acid" = non_genetic_sig_dha_diseases
)

library(ggvenn)
ggvenn(non_genetic_assoc, fill_color = c("blue4","#097969","#ff781f"),
       fill_alpha = 0.5,stroke_size = 1, set_name_size = 6, text_size = 4)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/NonGeneticAssociations_O3FA_O6FA_DHA_Protective_VennDiagram.pdf",height = 6,width = 6,dpi = 300,bg = "white")

#----------------------------------------------------------------------------------------

#--------------------Genetic association numbers----------------------
#--------------------------------O3FA-------------------------------------
# Number of total genetic associations O3FA
nrow(o3fa_prs_df[o3fa_prs_df$significance != "", ])

# Number of protective genetic associations O3FA
nrow(o3fa_prs_df[(o3fa_prs_df$significance != "") & (o3fa_prs_df$risk_prot == "Protective"), ])

# Unique disease groups for protective associations O3FA
unique(o3fa_prs_df[(o3fa_prs_df$significance != "") & (o3fa_prs_df$risk_prot == "Protective"), ]$category)

# Number of protective genetic associations O3FA
nrow(o3fa_prs_df[(o3fa_prs_df$significance != "") & (o3fa_prs_df$risk_prot == "Risk"), ])

# Unique disease groups for risk associations O3FA
unique(o3fa_prs_df[(o3fa_prs_df$significance != "") & (o3fa_prs_df$risk_prot == "Risk"), ]$category)

# Get significant O6FA associations
significant_prs_o3fa_df <- o3fa_prs_df[(o3fa_prs_df$significance != "") & (o3fa_prs_df$risk_prot == "Protective"), ]

st4_df <- significant_prs_o3fa_df[c("label","phenotype","category","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")]
colnames(st4_df) <- c("metabolite","disease","disease_group","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")
write.table(st4_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST4_genetic_O3FA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


# Make table of O3FA significant association category frequencies
category_significant_prs_o3fa_df <- count(significant_prs_o3fa_df,'category')
category_significant_prs_o3fa_df$from <- "Omega-3 fatty acids"

#--------------------------------O6FA-------------------------------------
# Number of total genetic associations O6FA
nrow(o6fa_prs_df[o6fa_prs_df$significance != "", ])

# Number of protective genetic associations O6FA
nrow(o6fa_prs_df[(o6fa_prs_df$significance != "") & (o6fa_prs_df$risk_prot == "Protective"), ])

# Unique disease groups for protective associations O6FA
unique(o6fa_prs_df[(o6fa_prs_df$significance != "") & (o6fa_prs_df$risk_prot == "Protective"), ]$category)

# Number of protective genetic associations O6FA
nrow(o6fa_prs_df[(o6fa_prs_df$significance != "") & (o6fa_prs_df$risk_prot == "Risk"), ])

# Unique disease groups for risk associations O6FA
unique(o6fa_prs_df[(o6fa_prs_df$significance != "") & (o6fa_prs_df$risk_prot == "Risk"), ]$category)

# Get significant O6FA associations
significant_prs_o6fa_df <- o6fa_prs_df[(o6fa_prs_df$significance != "") & (o6fa_prs_df$risk_prot == "Protective"), ]

st5_df <- significant_prs_o6fa_df[c("label","phenotype","category","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")]
colnames(st5_df) <- c("metabolite","disease","disease_group","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")
write.table(st5_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST5_genetic_O6FA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


# Make table of O6FA significant association category frequencies
category_significant_prs_o6fa_df <- count(significant_prs_o6fa_df,'category')
category_significant_prs_o6fa_df$from <- "Omega-6 fatty acids"

#--------------------------------DHA-------------------------------------
# Number of total genetic associations DHA
nrow(dha_prs_df[dha_prs_df$significance != "", ])

# Number of protective genetic associations DHA
nrow(dha_prs_df[(dha_prs_df$significance != "") & (dha_prs_df$risk_prot == "Protective"), ])

# Unique disease groups for protective associations DHA
unique(dha_prs_df[(dha_prs_df$significance != "") & (dha_prs_df$risk_prot == "Protective"), ]$category)

# Number of protective genetic associations DHA
nrow(dha_prs_df[(dha_prs_df$significance != "") & (dha_prs_df$risk_prot == "Risk"), ])

# Unique disease groups for risk associations DHA
unique(dha_prs_df[(dha_prs_df$significance != "") & (dha_prs_df$risk_prot == "Risk"), ]$category)

# Get significant DHA associations
significant_prs_dha_df <- dha_prs_df[(dha_prs_df$significance != "") & (dha_prs_df$risk_prot == "Protective"), ]

st6_df <- significant_prs_dha_df[c("label","phenotype","category","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")]
colnames(st6_df) <- c("metabolite","disease","disease_group","n_case_samples","n_control_samples","estimate","se","z","lo_ci","hi_ci","or","pval")
write.table(st6_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST6_genetic_DHA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


# Make table of O6FA significant association category frequencies
category_significant_prs_dha_df <- count(significant_prs_dha_df,'category')
category_significant_prs_dha_df$from <- "Docosahexaenoic acid"

# ------------------------Combine O3FA, O6FA, and DHA category counts for Circos plot----------------------------------
genetic_significant_protective_catcount_df <- rbind(category_significant_prs_dha_df,category_significant_prs_o6fa_df,category_significant_prs_o3fa_df)
genetic_significant_protective_catcount_df <- genetic_significant_protective_catcount_df[c('from','category','freq')]
colnames(genetic_significant_protective_catcount_df) <- c("from","to","value")

grid.col = c('Omega-3 fatty acids'='blue4','Omega-6 fatty acids'='#097969','Docosahexaenoic acid'='#ff781f',
             'Digestive'="#5A5A5A",'Dermatologic'='Grey','Circulatory system'='#5A5A5A','Endocrine/metabolic'="Grey",
             'Genitourinary'="#5A5A5A","Hematopoietic"="Grey","Infectious diseases"="#5A5A5A","Injuries & poisonings"="Grey",
             "Mental disorders"="#5A5A5A","Musculoskeletal"="Grey","Neoplasms"="#5A5A5A","Neurological"="Grey",
             "Respiratory"="#5A5A5A","Sense organs"="Grey","Symptoms"="#5A5A5A")

pdf("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/GeneticAssociations_O3FA_O6FA_DHA_Protective_ChordPlt.pdf",height = 20,width = 20)
chordDiagram(genetic_significant_protective_catcount_df,grid.col = grid.col,order=order)
dev.off()

#-------------------Venn Diagram for genetic associations-----------------------------
genetic_sig_o3fa_diseases <- as.list(significant_prs_o3fa_df$phenotype)
genetic_sig_o6fa_diseases <- as.list(significant_prs_o6fa_df$phenotype)
genetic_sig_dha_diseases <- as.list(significant_prs_dha_df$phenotype)

genetic_assoc <- list(
  "Omega-3 fatty acids" = genetic_sig_o3fa_diseases, 
  "Omega-6 fatty acids" = genetic_sig_o6fa_diseases, 
  "Docosahexaenoic acid" = genetic_sig_dha_diseases
)

library(ggvenn)
ggvenn(genetic_assoc, fill_color = c("blue4","#097969","#ff781f"),
       fill_alpha = 0.5,stroke_size = 1, set_name_size = 6, text_size = 4)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/GeneticAssociations_O3FA_O6FA_DHA_Protective_VennDiagram.pdf",height = 6,width = 6,dpi = 300,bg = "white")
#-------------------------------------------------------------------------------------

#-------------------Suggestive associations analysis-----------------------------
significant_o3fa_df <- significant_o3fa_df[c('label','phenotype','category','lo_ci','hi_ci','or','pval')]
colnames(significant_o3fa_df) <- c("label",'phenotype','category','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')
significant_o6fa_df <- significant_o6fa_df[c('label','phenotype','category','lo_ci','hi_ci','or','pval')]
colnames(significant_o6fa_df) <- c("label",'phenotype','category','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')
significant_dha_df <- significant_dha_df[c('label','phenotype','category','lo_ci','hi_ci','or','pval')]
colnames(significant_dha_df) <- c("label",'phenotype','category','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')

significant_prs_o3fa_df <- significant_prs_o3fa_df[c('label','phenotype','category','lo_ci','hi_ci','or','pval')]
colnames(significant_prs_o3fa_df) <- c("label",'phenotype','category','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')
significant_prs_o6fa_df <- significant_prs_o6fa_df[c('label','phenotype','category','lo_ci','hi_ci','or','pval')]
colnames(significant_prs_o6fa_df) <- c("label",'phenotype','category','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')
significant_prs_dha_df <- significant_prs_dha_df[c('label','phenotype','category','lo_ci','hi_ci','or','pval')]
colnames(significant_prs_dha_df) <- c("label",'phenotype','category','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')


significant_genetic_non_genetic_o3fa_df <- merge(significant_o3fa_df,significant_prs_o3fa_df,by="phenotype")
significant_genetic_non_genetic_o6fa_df <- merge(significant_o6fa_df,significant_prs_o6fa_df,by="phenotype")
significant_genetic_non_genetic_dha_df <- merge(significant_dha_df,significant_prs_dha_df,by="phenotype")

# Number of suggestive O3FA associations
nrow(significant_genetic_non_genetic_o3fa_df)
st7_df <- significant_genetic_non_genetic_o3fa_df[c("label.x","phenotype","category.x","non_genetic_lo_ci","non_genetic_hi_ci","non_genetic_or","non_genetic_pval","genetic_lo_ci","genetic_hi_ci","genetic_or","genetic_pval")]
colnames(st7_df) <- c("metabolite","disease","disease_group","non_genetic_lo_ci","non_genetic_hi_ci","non_genetic_or","non_genetic_pval","genetic_lo_ci","genetic_hi_ci","genetic_or","genetic_pval")
write.table(st7_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST7_suggestive_O3FA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

# Number of suggestive O6FA associations
nrow(significant_genetic_non_genetic_o6fa_df)
st8_df <- significant_genetic_non_genetic_o6fa_df[c("label.x","phenotype","category.x","non_genetic_lo_ci","non_genetic_hi_ci","non_genetic_or","non_genetic_pval","genetic_lo_ci","genetic_hi_ci","genetic_or","genetic_pval")]
colnames(st8_df) <- c("metabolite","disease","disease_group","non_genetic_lo_ci","non_genetic_hi_ci","non_genetic_or","non_genetic_pval","genetic_lo_ci","genetic_hi_ci","genetic_or","genetic_pval")
write.table(st8_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST8_suggestive_O6FA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

# Number of suggestive DHA associations
nrow(significant_genetic_non_genetic_dha_df)
st9_df <- significant_genetic_non_genetic_dha_df[c("label.x","phenotype","category.x","non_genetic_lo_ci","non_genetic_hi_ci","non_genetic_or","non_genetic_pval","genetic_lo_ci","genetic_hi_ci","genetic_or","genetic_pval")]
colnames(st9_df) <- c("metabolite","disease","disease_group","non_genetic_lo_ci","non_genetic_hi_ci","non_genetic_or","non_genetic_pval","genetic_lo_ci","genetic_hi_ci","genetic_or","genetic_pval")
write.table(st9_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST9_suggestive_DHA.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


category_significant_suggestive_o3fa_df <- count(significant_genetic_non_genetic_o3fa_df,'category.x')
category_significant_suggestive_o3fa_df$from <- "Omega-3 fatty acids"

category_significant_suggestive_o6fa_df <- count(significant_genetic_non_genetic_o6fa_df,'category.x')
category_significant_suggestive_o6fa_df$from <- "Omega-6 fatty acids"

category_significant_suggestive_dha_df <- count(significant_genetic_non_genetic_dha_df,'category.x')
category_significant_suggestive_dha_df$from <- "Docosahexaenoic acid"

# ------------------------Combine O3FA, O6FA, and DHA category counts for Circos plot----------------------------------
suggestive_significant_protective_catcount_df <- rbind(category_significant_suggestive_o3fa_df,category_significant_suggestive_o6fa_df,category_significant_suggestive_dha_df)
suggestive_significant_protective_catcount_df <- suggestive_significant_protective_catcount_df[c('from','category.x','freq')]
colnames(suggestive_significant_protective_catcount_df) <- c("from","to","value")

grid.col = c('Omega-3 fatty acids'='blue4','Omega-6 fatty acids'='#097969','Docosahexaenoic acid'='#ff781f',
             'Digestive'="#5A5A5A",'Dermatologic'='Grey','Circulatory system'='#5A5A5A','Endocrine/metabolic'="Grey",
             'Genitourinary'="#5A5A5A","Hematopoietic"="Grey","Infectious diseases"="#5A5A5A","Injuries & poisonings"="Grey",
             "Mental disorders"="#5A5A5A","Musculoskeletal"="Grey","Neoplasms"="#5A5A5A","Neurological"="Grey",
             "Respiratory"="#5A5A5A","Sense organs"="Grey","Symptoms"="#5A5A5A")

pdf("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/SuggestiveAssociations_O3FA_O6FA_DHA_Protective_ChordPlt.pdf",height = 20,width = 20)
chordDiagram(suggestive_significant_protective_catcount_df,grid.col = grid.col,order=order)
dev.off()

#-------------------Venn Diagram for suggestive associations-----------------------------
suggestive_sig_o3fa_diseases <- as.list(significant_genetic_non_genetic_o3fa_df$phenotype)
suggestive_sig_o6fa_diseases <- as.list(significant_genetic_non_genetic_o6fa_df$phenotype)
suggestive_sig_dha_diseases <- as.list(significant_genetic_non_genetic_dha_df$phenotype)

suggestive_assoc <- list(
  "Omega-3 fatty acids" = suggestive_sig_o3fa_diseases, 
  "Omega-6 fatty acids" = suggestive_sig_o6fa_diseases, 
  "Docosahexaenoic acid" = suggestive_sig_dha_diseases
)

ggvenn(suggestive_assoc, fill_color = c("blue4","#097969","#ff781f"),
       fill_alpha = 0.5,stroke_size = 1, set_name_size = 6, text_size = 4)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/SuggestiveAssociations_O3FA_O6FA_DHA_Protective_VennDiagram.pdf",height = 6,width = 6,dpi = 300,bg = "white")




# Create venn diagram for suggestive associations
# Create double bar plot for each disease group representing genetic and non-genetic protective p-values

#-----------------------------------------------------------
#-----------------------------------------------------------
colors  <- c("Circulatory system" = "blue4",
             "Congenital anomalies" = "royalblue1",
             "Dermatologic" = "blue4",
             "Digestive" = "royalblue1",
             "Endocrine/metabolic" = "blue4",
             "Genitourinary" = "royalblue1",
             "Hematopoietic" = "blue4",
             "Infectious diseases" = "royalblue1",
             "Injuries & poisonings" = "blue4",
             "Mental disorders" = "royalblue1",
             "Musculoskeletal" = "blue4",
             "Neoplasms" = "royalblue1",
             "Neurological" = "blue4",
             "Pregnancy complications" = "royalblue1",
             "Respiratory" = "blue4",
             "Sense organs" = "royalblue1",
             "Symptoms" = "blue4")

#------------------ Manhattan plot across disease groups for O3FA
ggplot(o3fa_df, aes(x=category,y=-log10(pval),color=category,label=significance,shape=risk_prot)) + coord_flip(ylim=c(0,15)) + geom_label_repel(aes(label=significance), segment.color="black", hjust=0,size=4, segment.size=0.25, nudge_x=0.5, direction="y") + 
  geom_point(aes(color = category,shape=risk_prot), size = 8,position = position_jitterdodge(dodge.width=0.01)) + 
  scale_color_manual(name="Cohort",values = colors) + geom_hline(yintercept=-log10(p_thresh), linetype="dashed", color = "grey", size=1) + 
  labs(x = "Disease Group", y="-log(p-value)",title = "Omega-3 fatty acids") + theme_classic() + scale_size_continuous(range = c(0.5,3)) +
  theme(axis.text.x = element_text(size = 20, hjust=1), axis.title = element_text(size = 20),legend.position="none") + 
  theme(text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/ManhattanPlot_Phecodes_MASTER_01_16_2023_Phecodes_PRS_O3FA_age_agesquare_sex_10PCS_glmresults.png",height = 15,width = 30,dpi = 300)


#------------------ Manhattan plot across disease groups for O6FA
colors2  <- c("Circulatory system" = "#097969",
              "Congenital anomalies" = "#AFE1AF",
              "Dermatologic" = "#097969",
              "Digestive" = "#AFE1AF",
              "Endocrine/metabolic" = "#097969",
              "Genitourinary" = "#AFE1AF",
              "Hematopoietic" = "#097969",
              "Infectious diseases" = "#AFE1AF",
              "Injuries & poisonings" = "#097969",
              "Mental disorders" = "#AFE1AF",
              "Musculoskeletal" = "#097969",
              "Neoplasms" = "#AFE1AF",
              "Neurological" = "#097969",
              "Pregnancy complications" = "#AFE1AF",
              "Respiratory" = "#097969",
              "Sense organs" = "#AFE1AF",
              "Symptoms" = "#097969")

ggplot(o6fa_df, aes(x=category,y=-log10(pval),color=category,label=significance,shape=risk_prot)) + coord_flip(ylim=c(0,15)) + geom_label_repel(aes(label=significance), segment.color="black", hjust=0,size=4, segment.size=0.25, nudge_x=0.5, direction="y") + 
  geom_point(aes(color = category,shape=risk_prot), size = 8,position = position_jitterdodge(dodge.width=0.01)) + 
  scale_color_manual(name="Cohort",values = colors2) + geom_hline(yintercept=-log10(p_thresh), linetype="dashed", color = "grey", size=1) + 
  labs(x = "Disease Group", y="-log(p-value)",title = "Omega-6 fatty acids") + theme_classic() + scale_size_continuous(range = c(0.5,3)) +
  theme(axis.text.x = element_text(size = 20, hjust=1), axis.title = element_text(size = 20),legend.position="none") +  
  theme(text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/ManhattanPlot_Phecodes_MASTER_01_16_2023_Phecodes_PRS_O6FA_age_agesquare_sex_10PCS_glmresults.png",height = 15,width = 30,dpi = 300)

#------------------ Manhattan plot across disease groups for DHA
colors3  <- c("Circulatory system" = "#ff781f",
              "Congenital anomalies" = "#ffaf7a",
              "Dermatologic" = "#ff781f",
              "Digestive" = "#ffaf7a",
              "Endocrine/metabolic" = "#ff781f",
              "Genitourinary" = "#ffaf7a",
              "Hematopoietic" = "#ff781f",
              "Infectious diseases" = "#ffaf7a",
              "Injuries & poisonings" = "#ff781f",
              "Mental disorders" = "#ffaf7a",
              "Musculoskeletal" = "#ff781f",
              "Neoplasms" = "#ffaf7a",
              "Neurological" = "#ff781f",
              "Pregnancy complications" = "#ffaf7a",
              "Respiratory" = "#ff781f",
              "Sense organs" = "#ffaf7a",
              "Symptoms" = "#ff781f")

ggplot(dha_df, aes(x=category,y=-log10(pval),color=category,label=significance,shape=risk_prot)) + coord_flip(ylim=c(0,15)) + geom_label_repel(aes(label=significance), segment.color="black", hjust=0,size=4, segment.size=0.25, nudge_x=0.5, direction="y") + 
  geom_point(aes(color = category,shape=risk_prot), size = 8,position = position_jitterdodge(dodge.width=0.01)) + 
  scale_color_manual(name="Cohort",values = colors3) + geom_hline(yintercept=-log10(p_thresh), linetype="dashed", color = "grey", size=1) + 
  labs(x = "Disease Group", y="-log(p-value)",title = "Docosahexaenoic Acid") + theme_classic() + scale_size_continuous(range = c(0.5,3)) +
  theme(axis.text.x = element_text(size = 20, hjust=1), axis.title = element_text(size = 20),legend.position="none") +  
  theme(text = element_text(family = "Helvetica",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20), panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/ManhattanPlot_Phecodes_MASTER_01_16_2023_Phecodes_PRS_DHA_age_agesquare_sex_10PCS_glmresults.png",height = 15,width = 30,dpi = 300)


#------------------ Comparison of O3FA and O6FA
o3fa_pruned_df1 <- o3fa_df[c("phenotype","pval")]
colnames(o3fa_pruned_df1) <- c("phenotype","o3fa_pval")

o6fa_pruned_df1 <- o6fa_df[c("phenotype","pval")]
colnames(o6fa_pruned_df1) <- c("phenotype","o6fa_pval")

dha_pruned_df1 <- dha_df[c("phenotype","pval")]
colnames(dha_pruned_df1) <- c("phenotype","dha_pval")

o3fa_o6fa_df <- merge(o3fa_pruned_df1,o6fa_pruned_df1,by="phenotype")
o3fa_o6fa_dha_df <- merge(o3fa_o6fa_df,dha_pruned_df1,by="phenotype")

o3fa_o6fa_dha_df <- o3fa_o6fa_dha_df[(o3fa_o6fa_dha_df$o3fa_pval < p_thresh) | (o3fa_o6fa_dha_df$o6fa_pval < p_thresh) | (o3fa_o6fa_dha_df$dha_pval < p_thresh), ]
either_or_df <- o3fa_o6fa_dha_df[c("phenotype")]

o3fa_pruned_df2 <- o3fa_df[c("phenotype","pval")]
o3fa_pruned_df2$label <- "Omega-3 fatty acids"

o6fa_pruned_df2 <- o6fa_df[c("phenotype","pval")]
o6fa_pruned_df2$label <- "Omega-6 fatty acids"

dha_pruned_df2 <- dha_df[c("phenotype","pval")]
dha_pruned_df2$label <- "Docosahexaenoic Acid"

o3fa_o6fa_dha_df2 <- rbind(o3fa_pruned_df2,o6fa_pruned_df2,dha_pruned_df2)

o3fa_o6fa_dha_df2 <- merge(either_or_df,o3fa_o6fa_dha_df2,by="phenotype")

o3fa_o6fa_dha_df3 <- merge(o3fa_o6fa_dha_df2,phecode_map_df,by="phenotype")

ggplot(o3fa_o6fa_dha_df3, aes(reorder(phenotype,(-log10(pval))), -log10(pval))) +   
  geom_bar(aes(fill = label), position = "dodge", stat="identity") + theme_classic() + 
  geom_hline(yintercept=-log10(p_thresh), linetype="dashed", color = "black", size=0.3) + 
  scale_fill_manual(values = c("#ff781f","blue4","#097969")) + 
  labs(x="Disease",y="-log(p-value)",fill="Group") + 
  theme(text = element_text(family = "Helvetica",size=20),legend.position="bottom") + facet_wrap(~category,scales="free") + coord_flip() 
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/PvalBarPlot_Phecodes_PRS_O3FA_VS_O6FA_VS_DHA_age_agesquare_sex_10PCS_glmresults.png",height = 50,width = 60,dpi = 300,limitsize = FALSE)

# Create venn diagrams 
sig_o3fa_df <- o3fa_df[o3fa_df$pval < p_thresh, ]
sig_o3fa_diseases <- as.list(sig_o3fa_df$phenotype)

sig_o6fa_df <- o6fa_df[o6fa_df$pval < p_thresh, ]
sig_o6fa_diseases <- as.list(sig_o6fa_df$phenotype)

sig_dha_df <- dha_df[dha_df$pval < p_thresh, ]
sig_dha_diseases <- as.list(sig_dha_df$phenotype)

x <- list(
  "Omega-3 fatty acids" = sig_o3fa_diseases, 
  "Omega-6 fatty acids" = sig_o6fa_diseases, 
  "Docosahexaenoic Acid" = sig_dha_diseases
)

library(ggvenn)
ggvenn(x, fill_color = c("blue4","#097969","#ff781f"),
       fill_alpha = 0.5,stroke_size = 1, set_name_size = 6, text_size = 4)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/VennDiagram_144Phecodes_PRS_O3FA_VS_O6FA_VS_DHA_age_agesquare_sex_10PCS_glmresults.png",height = 6,width = 6,dpi = 300,bg = "white")

