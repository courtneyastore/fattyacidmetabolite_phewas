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

def read_PRScs_ref(prs_cs_file):
    prs_cs_df = pd.read_csv(prs_cs_file,sep="\t")
    prs_cs_df['varID'] = prs_cs_df['CHR'].astype(str) + "-" + prs_cs_df['BP'].astype(str)
    prs_cs_df = prs_cs_df[['varID','SNP']]
    print(prs_cs_df)
    return prs_cs_df


# SNP          A1   A2   BETA      P
# A1 is the effect allele

def format_gwas(gwas_df):
    gwas_df['varID'] = gwas_df['chromosome'].astype(str) + "-" + gwas_df['base_pair_location'].astype(str)
    gwas_df = gwas_df[['varID','other_allele','effect_allele','beta','p_value']]
    gwas_df['other_allele'] = gwas_df['other_allele'].str.upper()
    gwas_df['effect_allele'] = gwas_df['effect_allele'].str.upper()

    gwas_df = gwas_df[['varID','effect_allele','other_allele','beta','p_value']]
    gwas_df.columns = ['varID','A1','A2','BETA','P']
    
    print(gwas_df)

    return gwas_df 

def merge_ref_gwas(gwas_df,prs_cs_df):
    merge_df = pd.merge(gwas_df,prs_cs_df,on="varID")
    merge_df = merge_df[['SNP','A1','A2','BETA','P']]
    print(merge_df)

    return merge_df

def write_gwas_file(final_df,final_file):
    
    final_df.to_csv(final_file,sep="\t",index=False,header=True)
    
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help = "Please input a GWAS file from Open GWAS: e.g., 28067908-GCST004131-EFO_0003767-build37.f.tsv")
    parser.add_argument("-p", help = "Please input PRScs reference SNP file: e.g., snpinfo_ukbb_hm3")
    args = parser.parse_args()

    GWASFile = read_gwas(args.f)
    PRScsFile = read_PRScs_ref(args.p)
    
    formatGWAS = format_gwas(GWASFile)

    mergeData = merge_ref_gwas(formatGWAS,PRScsFile)
    
    gwas_out_file = "format_PRScs_GWAS_"+ args.f
    
    writeGWASFile = write_gwas_file(mergeData,gwas_out_file)
    
if __name__ == "__main__":
    main()
