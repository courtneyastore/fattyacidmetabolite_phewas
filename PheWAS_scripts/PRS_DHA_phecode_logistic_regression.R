library(data.table)
library(dplyr)
library(stringr)

# ls ../processing_files/SplitPhecodesForParallelization/* | xargs -I {} -P 10 sh -c "Rscript metabolite_phecode_logistic_regression.R {}"

args <- commandArgs(trailingOnly = TRUE)

files_df <- as.data.frame(fread(args[1],sep="\t"))
files_lst <- as.list(files_df$phecode)

iter_ = basename(args[1])
iter_ = gsub("Cohort.","",iter_)
iter_ = gsub(".txt","",iter_)
out_file = "/storage/home/hcoda1/1/castore3/p-ggibson3-0/PUFA_PheWAS_Frontiers/result_files/11_25_2023_Phecodes_PRS_DHA_age_agesquare_sex_10PCS_glmresults_"
out_file = paste0(out_file,iter_,".tsv")

#files_lst <- list.files(path=args[1], pattern=".tsv", all.files=TRUE,full.names=TRUE)

master_df <- as.data.frame(fread("/storage/home/hcoda1/1/castore3/p-ggibson3-0/calculate_PRS/PRScs_scores/FinalPRSCalculations/met-d-dha_PRS.tsv",sep="\t"))
master_df$scalePRS <- scale(master_df$PRS)
master_df$eid <- as.character(master_df$ID)
head(master_df)


covar_df <- as.data.frame(fread("/storage/home/hcoda1/1/castore3/p-ggibson3-0/PUFA_PheWAS_Frontiers/processed_data/covariates_bodyweight_normalized.tsv",sep="\t"))
covar_df$eid <- paste(covar_df$eid,covar_df$eid,sep="_")
covar_df$eid <- as.character(covar_df$eid)


master_df <- merge(master_df,covar_df,by="eid")

#master_df <- master_df[c('eid','scalePRS','p31','p21022','p21022_square','p22009_a1','p22009_a2','p22009_a3','p22009_a4','p22009_a5','p22009_a6','p22009_a7','p22009_a8','p22009_a9','p22009_a10')]


#files_lst <- c('249','250')
glm_results_df <- NULL

for (phecode_id in files_lst){
	print(phecode_id)

	phecode_file = paste0("/storage/home/hcoda1/1/castore3/p-ggibson3-0/make_phecode_cohorts/scripts/PhecodeFiles/Cohort.",phecode_id,".tsv")

	case_control_cohort_df <- as.data.frame(fread(phecode_file,sep="\t"))

	# Rename column names.
	colnames(case_control_cohort_df) <- c('eid','disease_status')

	case_control_cohort_df$eid <- paste(case_control_cohort_df$eid,case_control_cohort_df$eid,sep="_")
	case_control_cohort_df$eid <- as.character(case_control_cohort_df$eid)

	# Merge independent variable data.
	case_control_cohort_df <- merge(master_df,case_control_cohort_df,by="eid",all.x=TRUE)
	#case_control_cohort_df <- case_control_cohort_df[complete.cases(case_control_cohort_df), ]
	
	# Remove excludes from cohort.
	case_control_cohort_df <- case_control_cohort_df[case_control_cohort_df$disease_status != "Excluded", ]

	case_control_cohort_df <- case_control_cohort_df[complete.cases(case_control_cohort_df$disease_status), ]
	
	# Recode case/control disease status to binary classification.
	case_control_cohort_df$disease_status <- ifelse(case_control_cohort_df$disease_status == "Case", 1, ifelse(case_control_cohort_df$disease_status == "Control", 0, NA))
	#case_control_cohort_df$disease_status <- ifelse(case_control_cohort_df$disease_status == "Case", 1,0)
	
	total_n_samples <- nrow(case_control_cohort_df)
    n_case_samples = nrow(case_control_cohort_df[case_control_cohort_df$disease_status == 1, ])
    n_control_samples = nrow(case_control_cohort_df[case_control_cohort_df$disease_status == 0, ])


	#case_control_cohort_df <- case_control_cohort_df[complete.cases(case_control_cohort_df), ]
	#print(nrow(case_control_cohort_df))
	
    if (n_case_samples < 50){
        print("Too few cases. Skipping.")
        next
    } else {
        glm_i <- glm(disease_status ~ scalePRS + p21022 + p21022_square + p31 + p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10, family = "binomial", data = case_control_cohort_df)

        estimate = coef(summary(glm_i))[2,1]
        se = coef(summary(glm_i))[2,2]
        z = coef(summary(glm_i))[2,3]
        or = exp(coef(glm_i))[2]
        pval = coef(summary(glm_i))[2,4]

        confint_ = exp(confint(glm_i))

        lo_ci = confint_[2,1]
        hi_ci = confint_[2,2]

        glm_results_df = rbind(glm_results_df, data.frame(phecode_id,total_n_samples,n_case_samples,n_control_samples,estimate,se,z,lo_ci,hi_ci,or,pval))

        #write.table(glm_results_df, file=out_file, quote=FALSE, sep='\t', col.names = TRUE,row.names=FALSE)
    }  
}

write.table(glm_results_df, file=out_file, quote=FALSE, sep='\t', col.names = TRUE,row.names=FALSE)

#phecode_df_indices = row.names(phecode_df)
#mclapply(phecode_df_indices, function(i) createPhecodeFiles(i), mc.cores = 32)

#write.table(glm_results_df, file='/storage/home/hcoda1/1/castore3/p-ggibson3-0/metabolitexenvironment_disease_project/result_files/Phecodes_Zscore_p23445_age_agesquare_sex_10PCS_glmresults.tsv', quote=FALSE, sep='\t', col.names = TRUE,row.names=FALSE)
