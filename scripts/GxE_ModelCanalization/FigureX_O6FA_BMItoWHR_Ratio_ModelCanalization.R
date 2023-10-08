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

phecode_map_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/phecode_definitions1.2.csv"))
phecode_map_df$phecode <- as.character(phecode_map_df$phecode)
phecode_map_df$category <- capFirst(phecode_map_df$category)

o3fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST8_suggestive_O3FA.tsv"))
o3fa_df <- o3fa_df[c('disease')]

o6fa_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST9_suggestive_O6FA.tsv"))
o6fa_df <- o6fa_df[c('disease')]

dha_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ST10_suggestive_DHA.tsv"))
dha_df <- dha_df[c('disease')]

all_suggestive_df <- rbind(o3fa_df,o6fa_df,dha_df)
all_suggestive_df <- unique(all_suggestive_df)

# Create delta departure df for O3FA BMI
o3fa_bmi_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_O6FAPGS_BMI_PrevPerc.tsv"))
o3fa_bmi_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_O6FAPGS_BMI_PrevPerc.tsv"))

o3fa_bmi_departure_df <- merge(o3fa_bmi_observed_df,o3fa_bmi_expected_df,by="phecode_id")
o3fa_bmi_departure_df <- o3fa_bmi_departure_df[c('phecode_id','delta','meanDelta','sdDelta','seDelta')]
o3fa_bmi_departure_df$phecode <- as.character(o3fa_bmi_departure_df$phecode_id)
o3fa_bmi_departure_df <- merge(o3fa_bmi_departure_df,phecode_map_df,by="phecode")
o3fa_bmi_departure_df <- o3fa_bmi_departure_df[o3fa_bmi_departure_df$category != "NULL", ]
o3fa_bmi_departure_df <- o3fa_bmi_departure_df[!(is.na(o3fa_bmi_departure_df$category)), ]
o3fa_bmi_departure_df <- o3fa_bmi_departure_df[o3fa_bmi_departure_df$category != "", ]
o3fa_bmi_departure_df <- o3fa_bmi_departure_df[c('phecode','delta','meanDelta','sdDelta','seDelta','phenotype','category')]
colnames(o3fa_bmi_departure_df) <- c('phecode','delta','meanDelta','sdDelta','seDelta','disease','category')
o3fa_bmi_departure_df <- merge(all_suggestive_df,o3fa_bmi_departure_df,by="disease")
o3fa_bmi_departure_df$departure <- (o3fa_bmi_departure_df$delta - o3fa_bmi_departure_df$meanDelta) / o3fa_bmi_departure_df$sdDelta # Departure = (deltaObs - deltaExp)/sd(deltaExp)

o3fa_bmi_departure_df <- o3fa_bmi_departure_df[c('disease','category','departure')]
colnames(o3fa_bmi_departure_df) <- c('disease','category','bmi_departure')

# Create delta departure df for O3FA WHR
o3fa_whr_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_O6FAPGS_WHR_PrevPerc.tsv"))
o3fa_whr_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_O6FAPGS_WHR_PrevPerc.tsv"))

o3fa_whr_departure_df <- merge(o3fa_whr_observed_df,o3fa_whr_expected_df,by="phecode_id")
o3fa_whr_departure_df <- o3fa_whr_departure_df[c('phecode_id','delta','meanDelta','sdDelta','seDelta')]
o3fa_whr_departure_df$phecode <- as.character(o3fa_whr_departure_df$phecode_id)
o3fa_whr_departure_df <- merge(o3fa_whr_departure_df,phecode_map_df,by="phecode")
o3fa_whr_departure_df <- o3fa_whr_departure_df[o3fa_whr_departure_df$category != "NULL", ]
o3fa_whr_departure_df <- o3fa_whr_departure_df[!(is.na(o3fa_whr_departure_df$category)), ]
o3fa_whr_departure_df <- o3fa_whr_departure_df[o3fa_whr_departure_df$category != "", ]
o3fa_whr_departure_df <- o3fa_whr_departure_df[c('phecode','delta','meanDelta','sdDelta','seDelta','phenotype','category')]
colnames(o3fa_whr_departure_df) <- c('phecode','delta','meanDelta','sdDelta','seDelta','disease','category')
o3fa_whr_departure_df <- merge(all_suggestive_df,o3fa_whr_departure_df,by="disease")
o3fa_whr_departure_df$departure <- (o3fa_whr_departure_df$delta - o3fa_whr_departure_df$meanDelta) / o3fa_whr_departure_df$sdDelta # Departure = (deltaObs - deltaExp)/sd(deltaExp)

o3fa_whr_departure_df <- o3fa_whr_departure_df[c('disease','departure')]
colnames(o3fa_whr_departure_df) <- c('disease','whr_departure')

merge_df <- merge(o3fa_whr_departure_df,o3fa_bmi_departure_df,by="disease")

merge_df$metabolite <- "O6FA"
merge_df <- merge_df[c('metabolite','disease','category','whr_departure','bmi_departure')]
#write.table(merge_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/tables/ForGreg_O6FA_WHR_BMI_departures.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names = TRUE)

nrow(merge_df[(merge_df$whr_departure > 2), ])
nrow(merge_df[(merge_df$whr_departure < -2), ])

nrow(merge_df[(merge_df$whr_departure > 1), ])
nrow(merge_df[(merge_df$whr_departure < -1), ])

merge_df$ratio <- merge_df$bmi_departure / merge_df$whr_departure


disease_colors  <- c("Circulatory system" = "#ef268f","Dermatologic" = "#d8946e","Digestive"="#55bb83", "Neurological" = "red",
                     "Endocrine/metabolic"="#2382de","Genitourinary"="#0ea9bb","Hematopoietic"="#b26f86","Infectious diseases"="#c8c2fd",
                     "Injuries & poisonings"="#6f0190","Mental disorders"="#bcd4b9","Musculoskeletal"="#a96e41","Neoplasms"="#a7b30f","Respiratory"="#ad4c45","Sense organs"="#dcac6d","Symptoms"="#598795")

merge_df <- merge_df[(merge_df$category == "Circulatory system") | (merge_df$category == "Digestive") | (merge_df$category == "Endocrine/metabolic"),]
#-----------------------------------------
ggplot(merge_df, aes(y=whr_departure,x=bmi_departure,color=category)) +
  scale_color_manual(name="Disease group",values = disease_colors) +
  geom_point(size=8) + theme_bw() + 
  labs(x = "BMI delta departure",y="WHR delta departure", title = "Omega-6 fatty acids") + 
  ylim(-6,4) + xlim(-6,4) + 
  geom_hline(yintercept=0, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=0, linetype="dashed", color = "grey") + 
  theme(axis.text.x = element_text(size = 20), text = element_text(family = "Helvetica",size=20),
        legend.title=element_text(size=20),legend.text=element_text(size=20),
        strip.text = element_text(size=20),legend.position="bottom",
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) +
  geom_text_repel(aes(label = ifelse(bmi_departure > 2 | bmi_departure < -2 | whr_departure > 2 | whr_departure < -2, disease, "")), 
                  size = 8, fontface = "bold")
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric_ModelCanalization_BMIWHRRatio/V3O6FA_BMI_DeltaDeparture_WHR_DeltaDeparture_Scatter_plt_Suggestive.pdf",height = 15,width = 15,dpi = 300,limitsize = FALSE)



o3fa_bmi_whr_departure_plt <- ggplot(merge_df, aes(reorder(disease,-ratio), ratio,color=category)) + 
  scale_color_manual(name="Disease group",values = disease_colors) +
  geom_point(size=5) + theme_bw() + geom_hline(yintercept=0) + coord_cartesian(ylim = c(-8, 8)) + 
  labs(x="Disease",y="BMI delta departure/WHR delta departure",title="Omega-6 fatty acids") +
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20))
o3fa_bmi_whr_departure_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric_ModelCanalization_BMIWHRRatio/O6FA_BMI_DeltaDeparture_WHR_DeltaDeparture_Suggestive.pdf",height = 15,width = 60,dpi = 300,limitsize = FALSE)

o3fa_bmi_whr_departure_plt <- ggplot(merge_df, aes(reorder(disease,-ratio), ratio,color=category)) + 
  scale_color_manual(name="Disease group",values = disease_colors) +
  geom_point(size=5) + theme_bw() + geom_hline(yintercept=0) + 
  labs(x="Disease",y="BMI delta departure/WHR delta departure",title="Omega-3 fatty acids") +
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20))
o3fa_bmi_whr_departure_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric_ModelCanalization_BMIWHRRatio/O6FA_BMI_DeltaDeparture_WHR_DeltaDeparture_Suggestive_all_points.pdf",height = 15,width = 60,dpi = 300,limitsize = FALSE)
