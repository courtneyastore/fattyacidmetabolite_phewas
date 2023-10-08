library(ggplot2)
library(forcats)
library(tidyverse)
library(dplyr)
library(data.table)
library(plyr)
library(stringr)

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
o3fa_bmi_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_O3FAPGS_BMI_PrevPerc.tsv"))
o3fa_bmi_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_O3FAPGS_BMI_PrevPerc.tsv"))

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

o3fa_bmi_departure_df <- o3fa_bmi_departure_df[c('disease','category','delta','departure')]
o3fa_bmi_departure_df$body_weight_measure <- "BMI"

# Create delta departure df for O3FA WHR
o3fa_whr_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_O3FAPGS_WHR_PrevPerc.tsv"))
o3fa_whr_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_O3FAPGS_WHR_PrevPerc.tsv"))

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

o3fa_whr_departure_df <- o3fa_whr_departure_df[c('disease','category','delta','departure')]
o3fa_whr_departure_df$body_weight_measure <- "WHR"

o3fa <- rbind(o3fa_whr_departure_df,o3fa_bmi_departure_df)
o3fa$metabolite <- "O3FA"

# Create delta departure df for O6FA BMI
o6fa_bmi_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_O6FAPGS_BMI_PrevPerc.tsv"))
o6fa_bmi_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_O6FAPGS_BMI_PrevPerc.tsv"))

o6fa_bmi_departure_df <- merge(o6fa_bmi_observed_df,o6fa_bmi_expected_df,by="phecode_id")
o6fa_bmi_departure_df <- o6fa_bmi_departure_df[c('phecode_id','delta','meanDelta','sdDelta','seDelta')]
o6fa_bmi_departure_df$phecode <- as.character(o6fa_bmi_departure_df$phecode_id)
o6fa_bmi_departure_df <- merge(o6fa_bmi_departure_df,phecode_map_df,by="phecode")
o6fa_bmi_departure_df <- o6fa_bmi_departure_df[o6fa_bmi_departure_df$category != "NULL", ]
o6fa_bmi_departure_df <- o6fa_bmi_departure_df[!(is.na(o6fa_bmi_departure_df$category)), ]
o6fa_bmi_departure_df <- o6fa_bmi_departure_df[o6fa_bmi_departure_df$category != "", ]
o6fa_bmi_departure_df <- o6fa_bmi_departure_df[c('phecode','delta','meanDelta','sdDelta','seDelta','phenotype','category')]
colnames(o6fa_bmi_departure_df) <- c('phecode','delta','meanDelta','sdDelta','seDelta','disease','category')
o6fa_bmi_departure_df <- merge(all_suggestive_df,o6fa_bmi_departure_df,by="disease")
o6fa_bmi_departure_df$departure <- (o6fa_bmi_departure_df$delta - o6fa_bmi_departure_df$meanDelta) / o6fa_bmi_departure_df$sdDelta # Departure = (deltaObs - deltaExp)/sd(deltaExp)

o6fa_bmi_departure_df <- o6fa_bmi_departure_df[c('disease','category','delta','departure')]
o6fa_bmi_departure_df$body_weight_measure <- "BMI"

# Create delta departure df for O6FA WHR
o6fa_whr_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_O6FAPGS_WHR_PrevPerc.tsv"))
o6fa_whr_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_O6FAPGS_WHR_PrevPerc.tsv"))

o6fa_whr_departure_df <- merge(o6fa_whr_observed_df,o6fa_whr_expected_df,by="phecode_id")
o6fa_whr_departure_df <- o6fa_whr_departure_df[c('phecode_id','delta','meanDelta','sdDelta','seDelta')]
o6fa_whr_departure_df$phecode <- as.character(o6fa_whr_departure_df$phecode_id)
o6fa_whr_departure_df <- merge(o6fa_whr_departure_df,phecode_map_df,by="phecode")
o6fa_whr_departure_df <- o6fa_whr_departure_df[o6fa_whr_departure_df$category != "NULL", ]
o6fa_whr_departure_df <- o6fa_whr_departure_df[!(is.na(o6fa_whr_departure_df$category)), ]
o6fa_whr_departure_df <- o6fa_whr_departure_df[o6fa_whr_departure_df$category != "", ]
o6fa_whr_departure_df <- o6fa_whr_departure_df[c('phecode','delta','meanDelta','sdDelta','seDelta','phenotype','category')]
colnames(o6fa_whr_departure_df) <- c('phecode','delta','meanDelta','sdDelta','seDelta','disease','category')
o6fa_whr_departure_df <- merge(all_suggestive_df,o6fa_whr_departure_df,by="disease")
o6fa_whr_departure_df$departure <- (o6fa_whr_departure_df$delta - o6fa_whr_departure_df$meanDelta) / o6fa_whr_departure_df$sdDelta # Departure = (deltaObs - deltaExp)/sd(deltaExp)

o6fa_whr_departure_df <- o6fa_whr_departure_df[c('disease','category','delta','departure')]
o6fa_whr_departure_df$body_weight_measure <- "WHR"

o6fa <- rbind(o6fa_bmi_departure_df,o6fa_whr_departure_df)
o6fa$metabolite <- "O6FA"

# Create delta departure df for DHA BMI
dha_bmi_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_DHAPGS_BMI_PrevPerc.tsv"))
dha_bmi_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_DHAPGS_BMI_PrevPerc.tsv"))

dha_bmi_departure_df <- merge(dha_bmi_observed_df,dha_bmi_expected_df,by="phecode_id")
dha_bmi_departure_df <- dha_bmi_departure_df[c('phecode_id','delta','meanDelta','sdDelta','seDelta')]
dha_bmi_departure_df$phecode <- as.character(dha_bmi_departure_df$phecode_id)
dha_bmi_departure_df <- merge(dha_bmi_departure_df,phecode_map_df,by="phecode")
dha_bmi_departure_df <- dha_bmi_departure_df[dha_bmi_departure_df$category != "NULL", ]
dha_bmi_departure_df <- dha_bmi_departure_df[!(is.na(dha_bmi_departure_df$category)), ]
dha_bmi_departure_df <- dha_bmi_departure_df[dha_bmi_departure_df$category != "", ]
dha_bmi_departure_df <- dha_bmi_departure_df[c('phecode','delta','meanDelta','sdDelta','seDelta','phenotype','category')]
colnames(dha_bmi_departure_df) <- c('phecode','delta','meanDelta','sdDelta','seDelta','disease','category')
dha_bmi_departure_df <- merge(all_suggestive_df,dha_bmi_departure_df,by="disease")
dha_bmi_departure_df$departure <- (dha_bmi_departure_df$delta - dha_bmi_departure_df$meanDelta) / dha_bmi_departure_df$sdDelta # Departure = (deltaObs - deltaExp)/sd(deltaExp)

dha_bmi_departure_df <- dha_bmi_departure_df[c('disease','category','delta','departure')]
dha_bmi_departure_df$body_weight_measure <- "BMI"

# Create delta departure df for DHA WHR
dha_whr_observed_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/ObservedDelta_DHAPGS_WHR_PrevPerc.tsv"))
dha_whr_expected_df <- as.data.frame(fread("/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/processing_files/ModelCanalization_BMIWHR/AvgSd_ExpectedDelta_DHAPGS_WHR_PrevPerc.tsv"))

dha_whr_departure_df <- merge(dha_whr_observed_df,dha_whr_expected_df,by="phecode_id")
dha_whr_departure_df <- dha_whr_departure_df[c('phecode_id','delta','meanDelta','sdDelta','seDelta')]
dha_whr_departure_df$phecode <- as.character(dha_whr_departure_df$phecode_id)
dha_whr_departure_df <- merge(dha_whr_departure_df,phecode_map_df,by="phecode")
dha_whr_departure_df <- dha_whr_departure_df[dha_whr_departure_df$category != "NULL", ]
dha_whr_departure_df <- dha_whr_departure_df[!(is.na(dha_whr_departure_df$category)), ]
dha_whr_departure_df <- dha_whr_departure_df[dha_whr_departure_df$category != "", ]
dha_whr_departure_df <- dha_whr_departure_df[c('phecode','delta','meanDelta','sdDelta','seDelta','phenotype','category')]
colnames(dha_whr_departure_df) <- c('phecode','delta','meanDelta','sdDelta','seDelta','disease','category')
dha_whr_departure_df <- merge(all_suggestive_df,dha_whr_departure_df,by="disease")
dha_whr_departure_df$departure <- (dha_whr_departure_df$delta - dha_whr_departure_df$meanDelta) / dha_whr_departure_df$sdDelta # Departure = (deltaObs - deltaExp)/sd(deltaExp)

dha_whr_departure_df <- dha_whr_departure_df[c('disease','category','delta','departure')]
dha_whr_departure_df$body_weight_measure <- "WHR"

dha <- rbind(dha_bmi_departure_df,dha_whr_departure_df)
dha$metabolite <- "DHA"

final_df <- rbind(dha,o3fa,o6fa)
final_df <- final_df[c('metabolite','disease','category','body_weight_measure','delta','departure')]

colnames(final_df) <- c('metabolite','disease','category','body_weight_measure','delta_observed','delta_departure')

final_df$delta_observed <- format(final_df$delta_observed, scientific = TRUE, digits = 2)
final_df$delta_departure <- format(final_df$delta_departure, scientific = TRUE, digits = 2)

write.table(final_df,file="/Users/courtneyastore/Dropbox (GaTech)/metabolitexenvironment_disease_project/PUFA_canalization_RshinyApp/PUFA_Canalization_table.tsv",quote=TRUE,sep="\t",row.names=FALSE,col.names = TRUE)








