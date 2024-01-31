#!/usr/bin/env python3

import argparse
import re
import os
import pandas as pd
import subprocess as sp

'''
DNAnexus
Path to WES plink files: /mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - interim 450k release
Files depicted as: ukb23149_c{1:22}_b*_v1.bim, ukb23149_c{1:22}_b*_v1.bed, ukb23149_c{1:22}_b*_v1.fam 
'''

# Read in Open GWAS summary statistic file: 28067908-GCST004131-EFO_0003767-build37.f.tsv
def read_gwas(gwas_file):
    gwas_df = pd.read_csv(gwas_file,sep="\t")
    print(gwas_df)
    return gwas_df

def format_gwas(gwas_df):
    gwas_df['varID'] = gwas_df['chromosome'].astype(str) + "-" + gwas_df['base_pair_location'].astype(str)
    gwas_df = gwas_df[['varID','other_allele','effect_allele','beta','p_value']]
    gwas_df['other_allele'] = gwas_df['other_allele'].str.upper()
    gwas_df['effect_allele'] = gwas_df['effect_allele'].str.upper()
    print(gwas_df)

    return gwas_df 

def write_gwas_file(pheno_df,pheno_file):
    
    pheno_df.to_csv(pheno_file,sep="\t",index=False,header=True)
    exit()
    proj_binary_cohorts_dir = "/mnt/project/cleaned_GWAS/"
    final_dir_file_dest = proj_binary_cohorts_dir + pheno_file

    # Make sure binary_cohorts directory exists in your project home directory. 
    if not os.path.exists(proj_binary_cohorts_dir):
        if not os.path.exists("cleaned_GWAS"):
            os.mkdir("cleaned_GWAS")
            sp.run(['dx','upload','-r','cleaned_GWAS'])
        else:
            sp.run(['dx','upload','-r','cleaned_GWAS'])
    else:
        print("cleaned_GWAS directory exists in your project home directory.")
        pass
    
    
    # Write your pheno file
    if not os.path.exists(final_dir_file_dest):   
        pheno_df.to_csv(pheno_file,sep="\t",index=False,header=True,na_rep='NA')
        sp.run(['dx','upload',pheno_file,'--path','cleaned_GWAS/'])
    else:
        print(final_dir_file_dest,"exists!")
        pass
    
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help = "Please input a GWAS file from Open GWAS: e.g., 28067908-GCST004131-EFO_0003767-build37.f.tsv")
    args = parser.parse_args()
    
    GWASFile = read_gwas(args.f)
    formatGWAS = format_gwas(GWASFile)
    
    gwas_out_file = "format_GWAS_"+ args.f
    
    writeGWASFile = write_gwas_file(formatGWAS,gwas_out_file)
    
if __name__ == "__main__":
    main()