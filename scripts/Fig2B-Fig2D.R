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
library(scales)
circos.par("track.height" = 0.4)
options(ggrepel.max.overlaps = Inf)

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

phecode_map_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/phecode_definitions1.2.csv"))
phecode_map_df$phecode <- as.character(phecode_map_df$phecode)
phecode_map_df$category <- capFirst(phecode_map_df$category)
phecode_category_count_df <- count(phecode_map_df,'category')

df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/result_files/MASTER_01_16_2023_Phecodes_PRS_O6FA_age_agesquare_sex_10PCS_glmresults.tsv"))
df$phecode <- as.character(df$phecode)
df <- merge(df,phecode_map_df,by="phecode")
df <- df[df$category != "NULL", ]
df <- df[!(is.na(df$category)), ]
category_count_df <- count(df,'category')
category_count_df <- category_count_df[category_count_df$category != "", ]

total_n_diseases = 1327

o3fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST7_suggestive_O3FA.tsv"))
o6fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST8_suggestive_O6FA.tsv"))
dha_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST9_suggestive_DHA.tsv"))

p_thresh
#-----------------------------------------------------------------------------------

# O3FA RR
n_sig_diseases = nrow(o3fa_df)

category_count_o3fa_df <- count(o3fa_df,'disease_group')

loop_lst <- as.list(category_count_o3fa_df$disease_group)

o3fa_p_thresh = 0.05/nrow(category_count_o3fa_df)

o3fa_final_df <- NULL
for (group in loop_lst){

  n_disease_group <- category_count_df[category_count_df$category == group, ][1,"freq"]
  n_sig_diseases_group <- category_count_o3fa_df[category_count_o3fa_df$disease_group == group, ][1,"freq"]
  
  n_sig_disease_n_sig_diseases_group = n_sig_diseases-n_sig_diseases_group
  total_n_diseases_n_diseases_group = total_n_diseases-n_disease_group
  
  cont_table = matrix(c(n_sig_diseases_group,n_sig_disease_n_sig_diseases_group,n_disease_group,total_n_diseases_n_diseases_group), nrow = 2)
  
  fishers_test = fisher.test(cont_table)
  
  rr = (n_sig_diseases_group/n_sig_diseases) / (n_disease_group/total_n_diseases)
  p = fishers_test$p.value
  o3fa_final_df = rbind(o3fa_final_df, data.frame(group,total_n_diseases,n_sig_diseases,n_disease_group,n_sig_diseases_group,rr,p))

}

o3fa_p_thresh_str_sig = paste0("p < ",as.character(scientific(o3fa_p_thresh, digits = 3)))
o3fa_p_thresh_str_sig

o3fa_final_df$Significance <- ifelse(o3fa_final_df$p<o3fa_p_thresh,"p < 3.57e-03","p > 3.57e-03")

colors  <- c("p < 3.57e-03" = "blue4",
             "p > 3.57e-03" = "royalblue1")

o3fa_plt <- ggplot(data=o3fa_final_df, aes(x=reorder(group,rr), y=rr,color=Significance)) + theme_bw() +
  scale_color_manual(name="Significance",values = colors) + 
  geom_point(size=8) + theme(legend.position="bottom",panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="RR",x="Disease group")
o3fa_plt

#ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/Fig1A_O3FA_Suggestive.pdf",height = 15,width = 80,dpi = 300,limitsize = FALSE)

#--------------------------------------------------------------------------------------
# O6FA RR
n_sig_diseases = nrow(o6fa_df)

category_count_o6fa_df <- count(o6fa_df,'disease_group')

o6fa_p_thresh = 0.05 / nrow(category_count_o6fa_df)

loop_lst <- as.list(category_count_o6fa_df$disease_group)

o6fa_final_df <- NULL
for (group in loop_lst){
  
  n_disease_group <- category_count_df[category_count_df$category == group, ][1,"freq"]
  n_sig_diseases_group <- category_count_o6fa_df[category_count_o6fa_df$disease_group == group, ][1,"freq"]
  
  n_sig_disease_n_sig_diseases_group = n_sig_diseases-n_sig_diseases_group
  total_n_diseases_n_diseases_group = total_n_diseases-n_disease_group
  
  cont_table = matrix(c(n_sig_diseases_group,n_sig_disease_n_sig_diseases_group,n_disease_group,total_n_diseases_n_diseases_group), nrow = 2)
  
  fishers_test = fisher.test(cont_table)
  
  rr = (n_sig_diseases_group/n_sig_diseases) / (n_disease_group/total_n_diseases)
  p = fishers_test$p.value
  o6fa_final_df = rbind(o6fa_final_df, data.frame(group,total_n_diseases,n_sig_diseases,n_disease_group,n_sig_diseases_group,rr,p))
  
}

#o6fa_p_thresh_str_sig = paste0("p < ",as.character(scientific(o6fa_p_thresh, digits = 3)))
#o6fa_p_thresh_str_notsig = paste0("p < ",as.character(scientific(o6fa_p_thresh, digits = 3)))

o6fa_final_df$Significance <- ifelse(o6fa_final_df$p<o6fa_p_thresh,"p < 4.17e-03","p > 4.17e-03")

colors  <- c("p < 4.17e-03" = "#097969",
             "p > 4.17e-03" = "#AFE1AF")

o6fa_plt <- ggplot(data=o6fa_final_df, aes(x=reorder(group,rr), y=rr,color=Significance)) + theme_bw() +
  scale_color_manual(name="Significance",values = colors) + 
  geom_point(size=8) + theme(legend.position="bottom",panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="RR",x="Disease group")
o6fa_plt


#--------------------------------------------------------------------------------------
# DHA RR
n_sig_diseases = nrow(dha_df)

category_count_dha_df <- count(dha_df,'disease_group')

dha_p_thresh = 0.05 / nrow(category_count_dha_df)

loop_lst <- as.list(category_count_dha_df$disease_group)

dha_final_df <- NULL
for (group in loop_lst){
  
  n_disease_group <- category_count_df[category_count_df$category == group, ][1,"freq"]
  n_sig_diseases_group <- category_count_dha_df[category_count_dha_df$disease_group == group, ][1,"freq"]
  
  n_sig_disease_n_sig_diseases_group = n_sig_diseases-n_sig_diseases_group
  total_n_diseases_n_diseases_group = total_n_diseases-n_disease_group
  
  cont_table = matrix(c(n_sig_diseases_group,n_sig_disease_n_sig_diseases_group,n_disease_group,total_n_diseases_n_diseases_group), nrow = 2)
  
  fishers_test = fisher.test(cont_table)
  
  rr = (n_sig_diseases_group/n_sig_diseases) / (n_disease_group/total_n_diseases)
  p = fishers_test$p.value
  dha_final_df = rbind(dha_final_df, data.frame(group,total_n_diseases,n_sig_diseases,n_disease_group,n_sig_diseases_group,rr,p))
  
}

dha_p_thresh_str_sig = paste0("p < ",as.character(scientific(dha_p_thresh, digits = 3)))
dha_p_thresh_str_sig

dha_final_df$Significance <- ifelse(dha_final_df$p<dha_p_thresh,"p < 3.33e-03","p > 3.33e-03")

colors  <- c("p < 3.33e-03" = "#ff781f",
             "p > 3.33e-03" = "#ffaf7a")

dha_plt <- ggplot(data=dha_final_df, aes(x=reorder(group,rr), y=rr,color=Significance)) + theme_bw() +
  scale_color_manual(name="Significance",values = colors) + 
  geom_point(size=8) + theme(legend.position="bottom",panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) + 
  labs(y="RR",x="Disease group")
dha_plt

final_plt <- grid.arrange(o3fa_plt, o6fa_plt, dha_plt, nrow = 1)
final_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/DiseaseGroup_RR.pdf",final_plt,height = 7,width = 22,dpi = 300,limitsize = FALSE)

