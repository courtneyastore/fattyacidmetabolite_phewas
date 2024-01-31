#!/usr/bin/env python3

import argparse
import re
import os
import numpy as np
import pandas as pd
import subprocess as sp

'''
DNAnexus
Path to WES plink files: /mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - interim 450k release
Files depicted as: ukb23149_c{1:22}_b*_v1.bim, ukb23149_c{1:22}_b*_v1.bed, ukb23149_c{1:22}_b*_v1.fam 
'''

# Read in Open GWAS summary statistic file: 28067908-GCST004131-EFO_0003767-build37.f.tsv
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	met-d-Omega_3
def read_gwas(f):
    all_data = pd.read_csv(f,sep='\t',low_memory=False,comment="#",names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GWAS"])
    df = pd.DataFrame(all_data)

    current_header = df.columns.values
    gwas_info = df['FORMAT'].iloc[0].split(":")

    df[gwas_info] = df['GWAS'].str.split(':', n=len(gwas_info), expand=True)

    df = df[['ID','REF','ALT','ES','LP']]

    df['P'] = df['LP'].astype(float) * -1
    df['P'] = np.power(10,df['P'])

    df = df[['ID','REF','ALT','ES','P']]

    df.columns = ['SNP','A2','A1','BETA','P']

    print(df)
    return df

def read_PRScs_ref(prs_cs_file):
    prs_cs_df = pd.read_csv(prs_cs_file,sep="\t")
    prs_cs_df['varID'] = prs_cs_df['CHR'].astype(str) + "-" + prs_cs_df['BP'].astype(str)
    prs_cs_df = prs_cs_df[['varID','SNP']]
    print(prs_cs_df)

    return prs_cs_df


# SNP          A1   A2   BETA      P
# A1 is the effect allele


def merge_ref_gwas(gwas_df,prs_cs_df):
    merge_df = pd.merge(gwas_df,prs_cs_df,on="SNP")
    merge_df = merge_df[['SNP','A1','A2','BETA','P']]
    merge_df = merge_df.drop_duplicates()
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

    mergeData = merge_ref_gwas(GWASFile,PRScsFile)
    
    gwas_out_file = "format_PRScs_GWAS_"+ (args.f).replace(".gz","")
    
    writeGWASFile = write_gwas_file(mergeData,gwas_out_file)
    
if __name__ == "__main__":
    main()
