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

o3fa_bmi_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_O6FAPGS_BMI_PrevPerc.tsv"))
o3fa_bmi_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_O6FAPGS_BMI_PrevPerc.tsv"))

# Create delta departure df for O3FA BMI
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


# Comparing observed vs expected delta for O3FA BMI 
o3fa_bmi_observed_df$label <- "Observed"
o3fa_bmi_observed_df$phecode <- as.character(o3fa_bmi_observed_df$phecode_id)
o3fa_bmi_observed_df <- o3fa_bmi_observed_df[c('phecode','delta','label')]
colnames(o3fa_bmi_observed_df) <- c('phecode','delta','label')

o3fa_bmi_expected_df$label <- "Expected"
o3fa_bmi_expected_df$phecode <- as.character(o3fa_bmi_expected_df$phecode_id)
o3fa_bmi_expected_df <- o3fa_bmi_expected_df[c('phecode','meanDelta','label')]
colnames(o3fa_bmi_expected_df) <- c('phecode','delta','label')

o3fa_bmi_df <- rbind(o3fa_bmi_observed_df,o3fa_bmi_expected_df)
o3fa_bmi_df <- merge(o3fa_bmi_df,phecode_map_df,by="phecode")
o3fa_bmi_df <- o3fa_bmi_df[o3fa_bmi_df$category != "NULL", ]
o3fa_bmi_df <- o3fa_bmi_df[!(is.na(o3fa_bmi_df$category)), ]
o3fa_bmi_df <- o3fa_bmi_df[o3fa_bmi_df$category != "", ]
o3fa_bmi_df <- o3fa_bmi_df[c('phecode','phenotype','category','delta','label')]
colnames(o3fa_bmi_df) <- c('phecode','disease','category','delta','label')

o3fa_bmi_df <- merge(all_suggestive_df,o3fa_bmi_df,by="disease")

colors  <- c("Observed" = "#097969",
             "Expected" = "#AFE1AF")

o3fa_whr_plt <- ggplot(o3fa_bmi_df, aes(reorder(disease,-delta), delta, fill = label)) + 
  scale_fill_manual(name="Delta type",values = colors) +
  geom_bar(stat="identity", position = "dodge") + theme_bw() + 
  labs(x="Disease",y="Delta",title="E. Omega-6 fatty acids and BMI") +
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20))
o3fa_whr_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric_ModelCanalization/V2Fig7E_O6FA_BMI_ObservedExpectedDelta_Suggestive.pdf",height = 15,width = 60,dpi = 300,limitsize = FALSE)

disease_colors  <- c("Circulatory system" = "#ef268f","Dermatologic" = "#d8946e","Digestive"="#55bb83", "Neurological" = "red",
                     "Endocrine/metabolic"="#2382de","Genitourinary"="#0ea9bb","Hematopoietic"="#b26f86","Infectious diseases"="#c8c2fd",
                     "Injuries & poisonings"="#6f0190","Mental disorders"="#bcd4b9","Musculoskeletal"="#a96e41","Neoplasms"="#a7b30f","Respiratory"="#ad4c45","Sense organs"="#dcac6d","Symptoms"="#598795")

o3fa_whr_departure_plt <- ggplot(o3fa_bmi_departure_df, aes(reorder(disease,-departure), departure,color=category)) + 
  scale_color_manual(name="Disease group",values = disease_colors) +
  geom_point(size=5) + theme_bw() + geom_hline(yintercept=0) + 
  labs(x="Disease",y="Delta departure",title="F. Omega-6 fatty acids and BMI (Departure)") +
  theme(axis.text.x = element_text(size = 20,angle = 45, hjust=1), text = element_text(family = "Helvetica",size=20),legend.title=element_text(size=20),legend.text=element_text(size=20),strip.text = element_text(size=20),legend.position="bottom",axis.text.y = element_text(size = 20), axis.title = element_text(size = 20))
o3fa_whr_departure_plt
ggsave("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/figures/MetaboliteXAnthropometric_ModelCanalization/V2Fig7F_O6FA_BMI_DeltaDeparture_Suggestive.pdf",height = 15,width = 60,dpi = 300,limitsize = FALSE)
