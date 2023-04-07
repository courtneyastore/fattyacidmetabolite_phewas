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

o3fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST7_suggestive_O3FA.tsv"))
o3fa_suggestive_disease_df <- o3fa_df[c('disease')]
o6fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST8_suggestive_O6FA.tsv"))
o6fa_suggestive_disease_df <- o6fa_df[c('disease')]
dha_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST9_suggestive_DHA.tsv"))
dha_suggestive_disease_df <- dha_df[c('disease')]
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
colnames(mr_o3fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
mr_o3fa_df$label <- "Genetic (Mendelian randomization)"
mr_o3fa_df <- mr_o3fa_df[(mr_o3fa_df$p < 0.05) & (mr_o3fa_df$or < 1), ]
mr_o3fa_df <- merge(mr_o3fa_df,o3fa_suggestive_disease_df,by="disease")

# Number of O3FA MR disease associations
nrow(mr_o3fa_df)

o3fa_disease_df <- mr_o3fa_df[c('disease')]


# Nongenetic
nongenetic_o3fa_df <- o3fa_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_o3fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_o3fa_df$label <- "Nongenetic"

# Genetic
genetic_o3fa_df <- o3fa_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_o3fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_o3fa_df$label <- "Genetic (PRS)"

final_o3fa_df <- rbind(nongenetic_o3fa_df,genetic_o3fa_df,mr_o3fa_df)

final_o3fa_df <- merge(o3fa_disease_df,final_o3fa_df,by="disease")

colors  <- c("Nongenetic" = "blue4",
             "Genetic (PRS)" = "royalblue1",
             "Genetic (Mendelian randomization)" = "#ADD8E6")

# Genetic vs non genetic vs MR
ggplot(data=final_o3fa_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=2,position = position_dodge(width = 0.3)) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=5) + facet_grid(. ~ disease_group, scales = "free_x", space='free') + ylim(0.75,1)
#ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig3A_O3FA_Suggestive_MR.pdf",height = 8,width = 16,dpi = 300,limitsize = FALSE)



#----------------------------------------------

# MR
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
colnames(mr_o6fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
mr_o6fa_df$label <- "Genetic (Mendelian randomization)"
mr_o6fa_df <- mr_o6fa_df[(mr_o6fa_df$p < 0.05) & (mr_o6fa_df$or < 1), ]
mr_o6fa_df <- merge(mr_o6fa_df,o6fa_suggestive_disease_df,by="disease")

# Number of O3FA MR disease associations
nrow(mr_o6fa_df)

o6fa_disease_df <- mr_o6fa_df[c('disease')]


# Nongenetic
nongenetic_o6fa_df <- o6fa_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_o6fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_o6fa_df$label <- "Nongenetic"

# Genetic
genetic_o6fa_df <- o6fa_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_o6fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_o6fa_df$label <- "Genetic (PRS)"

final_o6fa_df <- rbind(nongenetic_o6fa_df,genetic_o6fa_df,mr_o6fa_df)

final_o6fa_df <- merge(o6fa_disease_df,final_o6fa_df,by="disease")

colors  <- c("Nongenetic" = "#097969",
             "Genetic (PRS)" = "#AFE1AF",
             "Genetic (Mendelian randomization)"="#50C878")

# Genetic vs non genetic vs MR
ggplot(data=final_o6fa_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=2,position = position_dodge(width = 0.5)) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=5) + facet_grid(. ~ disease_group, scales = "free_x", space='free') + ylim(0.1,1)

#ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig3A_O6FA_Suggestive_MR.pdf",height = 12,width = 24,dpi = 300,limitsize = FALSE)


#----------------------------------
# MR
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
colnames(mr_dha_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
mr_dha_df$label <- "Genetic (Mendelian randomization)"
mr_dha_df <- mr_dha_df[(mr_dha_df$p < 0.05) & (mr_dha_df$or < 1), ]
mr_dha_df <- merge(mr_dha_df,dha_suggestive_disease_df,by="disease")

# Number of DHA MR disease associations
nrow(mr_dha_df)

dha_disease_df <- mr_dha_df[c('disease')]


# Nongenetic
nongenetic_dha_df <- dha_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_dha_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_dha_df$label <- "Nongenetic"

# Genetic
genetic_dha_df <- dha_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_dha_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_dha_df$label <- "Genetic (PRS)"

final_dha_df <- rbind(nongenetic_dha_df,genetic_dha_df,mr_dha_df)

final_dha_df <- merge(dha_disease_df,final_dha_df,by="disease")

colors  <- c("Nongenetic" = "#ff781f",
             "Genetic (PRS)" = "#ffaf7a",
             "Genetic (Mendelian randomization)"="#F88379")

# Genetic vs non genetic vs MR
ggplot(data=final_dha_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=2,position = position_dodge(width = 0.5)) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=5) + facet_grid(. ~ disease_group, scales = "free_x", space='free') + ylim(0.5,1)

#ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig3A_DHA_Suggestive_MR.pdf",height = 12,width = 37,dpi = 300,limitsize = FALSE)

#-----------------------------------------------------------------------
# Final diseases for PRS calculation.
gwas_id_map_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/All_SuggestiveDiseases.txt"))
final_mr_diseases_df <- rbind(o3fa_disease_df,o6fa_disease_df,dha_disease_df)
final_mr_diseases_df <- unique(final_mr_diseases_df)
final_mr_diseases_df <- merge(gwas_id_map_df,final_mr_diseases_df,by="disease")
write.table(final_mr_diseases_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/All_SuggestiveDiseases_MRfindings_forPRS.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
#-----------------------------------------------------------------------










# O6FA
nongenetic_o6fa_df <- o6fa_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_o6fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_o6fa_df$label <- "Nongenetic"

genetic_o6fa_df <- o6fa_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_o6fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_o6fa_df$label <- "Genetic"

final_o6fa_df <- rbind(nongenetic_o6fa_df,genetic_o6fa_df)

colors  <- c("Nongenetic" = "#097969",
             "Genetic" = "#AFE1AF")

# OR genetic and non-genetic associations O3FA
ggplot(data=final_o6fa_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=1) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=5) + facet_grid(. ~ disease_group, scales = "free_x", space='free') + ylim(0.2,1)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig2A_O6FA_Suggestive.pdf",height = 15,width = 80,dpi = 300,limitsize = FALSE)

# Genetic vs non genetic -log10 p O3FA
ggplot(data=o6fa_df,aes(x=-log10(non_genetic_pval),y=-log10(non_genetic_pval),color=disease_group,label=disease)) + 
  labs(x="-log(non-genetic p-value)",y="-log(genetic p-value)") + theme_bw() + 
  theme(text = element_text(family = "Helvetica",size=20),strip.text = element_text(size=20),legend.position="none",axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) +
  geom_point(size=2) + geom_text_repel() + facet_wrap(~disease_group,scales="free")
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig2B_O6FA_Suggestive.pdf",height = 20,width = 20,dpi = 300,limitsize = FALSE)

# DHA
nongenetic_dha_df <- dha_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_dha_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_dha_df$label <- "Nongenetic"

genetic_dha_df <- dha_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_dha_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_dha_df$label <- "Genetic"

final_dha_df <- rbind(nongenetic_dha_df,genetic_dha_df)

colors  <- c("Nongenetic" = "#ff781f",
             "Genetic" = "#ffaf7a")

# OR genetic and non-genetic associations O3FA
ggplot(data=final_dha_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=1) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=5) + facet_grid(. ~ disease_group, scales = "free_x", space='free') + ylim(0.3,1)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig3A_DHA_Suggestive.pdf",height = 15,width = 80,dpi = 300,limitsize = FALSE)







# Genetic vs non genetic -log10 p O3FA
ggplot(data=dha_df,aes(x=-log10(non_genetic_pval),y=-log10(non_genetic_pval),color=disease_group,label=disease)) + 
  labs(x="-log(non-genetic p-value)",y="-log(genetic p-value)") + theme_bw() + 
  theme(text = element_text(family = "Helvetica",size=20),strip.text = element_text(size=20),legend.position="none",axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) +
  geom_point(size=2) + geom_text_repel() + facet_wrap(~disease_group,scales="free")
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig3B_DHA_Suggestive.pdf",height = 20,width = 20,dpi = 300,limitsize = FALSE)







#--------------------------Archive-------------------------------
# OR genetic and non-genetic associations O3FA
ggplot(data=final_o3fa_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=2) + theme(text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=2) + coord_flip() + facet_wrap(~disease_group,scales="free") + ylim(0.4,1)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1A_O3FA_Suggestive.pdf",height = 30,width = 50,dpi = 300,limitsize = FALSE)

# OR genetic and non-genetic associations O3FA
ggplot(data=final_o3fa_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=2) + theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=2) + facet_wrap(~disease_group,nrow=1,scales="free_x") + ylim(0.4,1)
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1A_O3FA_Suggestive.pdf",height = 20,width = 80,dpi = 300,limitsize = FALSE)

