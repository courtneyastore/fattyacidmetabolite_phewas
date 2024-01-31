Step 1. Download VCF file from Open GWAS
wget https://gwas.mrcieu.ac.uk/files/met-d-Omega_3/met-d-Omega_3.vcf.gz

Step 2. Format GWAS summary statistics from Open GWAS 
cd /storage/home/hcoda1/1/castore3/p-ggibson3-0/calculate_PRS/raw_GWAS
python3 ../scripts/format_OpenGWAS_GRCh37_sumstats_fpr_PRScs.py -f met-d-Omega_3.vcf.gz -p ~/scratch/ldblk_ukbb_eur/snpinfo_ukbb_hm3
mv format_PRScs_GWAS_met-d-Omega_3.vcf /storage/home/hcoda1/1/castore3/p-ggibson3-0/calculate_PRS/format_PRScs/

Step 3. Run PRS-cs
echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 | tr " " "\n" | xargs -I {} -P 11 sh -c "python ../../../PRScs/PRScs.py --ref_dir=/storage/home/hcoda1/1/castore3/scratch/ldblk_ukbb_eur --bim_prefix=/storage/home/hcoda1/1/castore3/scratch/UKBB/GENOTYPE/BED_AFTER_QC/UKB_BGENQC_TO_BED/UKB_BED_AFTER_QC_chr{} --sst_file=/storage/home/hcoda1/1/castore3/p-ggibson3-0/calculate_PRS/format_PRScs/format_PRScs_GWAS_met-d-Omega_3.vcf --chrom={} --phi=1e-2 --out_dir=Chr{} --n_gwas=114999"

Step 4. Calculate PRS using plink2
/gpfs/pace1/project/bio-gibson/snagpal3/Tools/plink2_Mar19 --bfile /gpfs/scratch1/3/snagpal3/UKB_TRAITS_PRUNING_PRS/EDUC_ATTAINMENT_Prune_PRS/MERGED/MERGED_EDUC --keep "$WB_ids" --score /gpfs/scratch1/3/snagpal3/UKB_TRAITS_PRUNING_PRS/EDUC_ATTAINMENT_Prune_PRS/PRS/Profile_5e08_LD1.txt cols=+scoresums list-variants --out /gpfs/scratch1/3/snagpal3/UKB_TRAITS_PRUNING_PRS/EDUC_ATTAINMENT_Prune_PRS/PRS/PRS_5e-08

echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 | tr " " "\n" | xargs -I {} -P 11 plink2 --bfile ~/scratch/UKBB/GENOTYPE/BED_AFTER_QC/UKB_BGENQC_TO_BED/UKB_BED_AFTER_QC_chr{} --score Chr{}_pst_eff_a1_b0.5_phi1e-02_chr{}.txt 2 4 6 cols=+scoresums list-variants --out Chr{}_UKBB_PRScs_deLange

Step 5. Get the final PRS score across all chromosomes
python3 calculate_PRS_across_chromosomes.py -d ../PRScs_scores/deLange_IBD/ -o ../PRScs_scores/FinalPRSCalculations/
