library(dplyr)

# get a list of all file names in the directory
files <- list.files(path = "/storage/home/hcoda1/1/castore3/p-ggibson3-0/calculate_PRS/PRScs_scores/FinalPRSCalculations/", pattern = ".tsv", full.names = TRUE)

# read the first file as the base dataset
base_df <- read.csv(files[1],sep="\t")
base_df$PRS <- scale(base_df$PRS)
file_name <- gsub("\\..*", "", basename(files[1]))
colnames(base_df) <- ifelse(colnames(base_df) == "ID", colnames(base_df), paste(colnames(base_df), file_name, sep = "_"))

# loop over the remaining files and merge them with the base dataset
for (i in 2:length(files)) {
  file <- read.csv(files[i],sep="\t")
  file$PRS <- scale(file$PRS)
  print(head(file))
  # get the name of the file without the file extension
  file_name <- gsub("\\..*", "", basename(files[i]))
  
  # rename the columns in the file
  colnames(file) <- ifelse(colnames(file) == "ID", colnames(file), paste(colnames(file), file_name, sep = "_"))

  base_df <- full_join(base_df, file, by = "ID")

}

# write the merged dataset to a new file
write.table(base_df, "/storage/home/hcoda1/1/castore3/p-ggibson3-0/calculate_PRS/PRScs_scores/V2ALLmergedPRS_UKBB.tsv",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
