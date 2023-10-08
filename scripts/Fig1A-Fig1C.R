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


o3fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST8_suggestive_O3FA.tsv"))
o6fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST9_suggestive_O6FA.tsv"))
dha_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST10_suggestive_DHA.tsv"))

# O3FA
nongenetic_o3fa_df <- o3fa_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_o3fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_o3fa_df$label <- "Nongenetic"

genetic_o3fa_df <- o3fa_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_o3fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_o3fa_df$label <- "Genetic (PGS)"

final_o3fa_df <- rbind(nongenetic_o3fa_df,genetic_o3fa_df)

colors  <- c("Nongenetic" = "blue4",
             "Genetic (PGS)" = "royalblue1")


# OR genetic and non-genetic associations O3FA
ggplot(data=final_o3fa_df, aes(y=disease, x=or,xmin=lo_ci, xmax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=1) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_vline(xintercept=1, lty=5) + facet_grid(disease_group ~ ., scales = "free_y", space='free') + xlim(0.4,1) 
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1A_O3FA_Suggestive.pdf",height = 30,width = 15,dpi = 300,limitsize = FALSE)



# O6FA
nongenetic_o6fa_df <- o6fa_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_o6fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_o6fa_df$label <- "Nongenetic"

genetic_o6fa_df <- o6fa_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_o6fa_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_o6fa_df$label <- "Genetic (PGS)"

final_o6fa_df <- rbind(nongenetic_o6fa_df,genetic_o6fa_df)

colors  <- c("Nongenetic" = "#097969",
             "Genetic (PGS)" = "#AFE1AF")

# OR genetic and non-genetic associations O3FA
ggplot(data=final_o6fa_df, aes(y=disease, x=or,xmin=lo_ci, xmax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=1) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_vline(xintercept=1, lty=5) + facet_grid(disease_group ~ ., scales = "free_y", space='free') + xlim(0.2,1) 
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1B_O6FA_Suggestive.pdf",height = 40,width = 15,dpi = 300,limitsize = FALSE)



# DHA
nongenetic_dha_df <- dha_df[c('metabolite','disease','disease_group','non_genetic_lo_ci','non_genetic_hi_ci','non_genetic_or','non_genetic_pval')]
colnames(nongenetic_dha_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
nongenetic_dha_df$label <- "Nongenetic"

genetic_dha_df <- dha_df[c('metabolite','disease','disease_group','genetic_lo_ci','genetic_hi_ci','genetic_or','genetic_pval')]
colnames(genetic_dha_df) <- c('metabolite','disease','disease_group','lo_ci','hi_ci','or','p')
genetic_dha_df$label <- "Genetic (PGS)"

final_dha_df <- rbind(nongenetic_dha_df,genetic_dha_df)

colors  <- c("Nongenetic" = "#ff781f",
             "Genetic (PGS)" = "#ffaf7a")

# OR genetic and non-genetic associations O3FA
ggplot(data=final_o6fa_df, aes(y=disease, x=or,xmin=lo_ci, xmax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=1) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_vline(xintercept=1, lty=5) + facet_grid(disease_group ~ ., scales = "free_y", space='free') + xlim(0.2,1) 
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1C_DHA_Suggestive.pdf",height = 40,width = 15,dpi = 300,limitsize = FALSE)













#--------------------------Archive-------------------------------
# OR genetic and non-genetic associations O3FA
ggplot(data=final_o3fa_df, aes(x=disease, y=or,ymin=lo_ci, ymax=hi_ci,color=label)) + theme_bw() + scale_color_manual(name="Association type",values = colors) + 
  geom_pointrange(size=1) + theme(strip.background.x = element_rect(fill = "grey", linetype = "solid",color = "black", linewidth = 0), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="OR (95% CI)",x="",color="Association type") + geom_hline(yintercept=1, lty=5) + facet_grid(. ~ disease_group, scales = "free_x", space='free') + ylim(0.4,1) 
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1A_O3FA_Suggestive.pdf",height = 15,width = 80,dpi = 300,limitsize = FALSE)


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

# Genetic vs non genetic -log10 p O3FA
ggplot(data=o3fa_df,aes(x=-log10(non_genetic_pval),y=-log10(non_genetic_pval),color=disease_group,label=disease)) + 
  labs(x="-log(non-genetic p-value)",y="-log(genetic p-value)") + theme_bw() + 
  theme(text = element_text(family = "Helvetica",size=20),strip.text = element_text(size=20),legend.position="none",axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) +
  geom_point(size=2) + geom_text_repel() + facet_wrap(~disease_group,scales="free")
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1B_O3FA_Suggestive.pdf",height = 20,width = 20,dpi = 300,limitsize = FALSE)
