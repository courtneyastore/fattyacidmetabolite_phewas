#!/usr/bin/env python3

import argparse
import re
import os
import pandas as pd
import subprocess as sp

# Get all .sscore files in directory.
def score_dir(d):
    score_lst = [i for i in os.listdir(d) if ".sscore" in i and ".sscore.vars" not in i]

    return score_lst

# Loop over all of the plink score files and merge to one file per chr.
def compute_PRS(d,score_lst,o):
    # Define output file name.
    out_file_name = str(o) + os.path.basename(os.path.normpath(d)) + "_PRS.tsv"
    
    # Create NULL df.
    merged_df = pd.DataFrame()
    for i in score_lst:
        i_dir = str(d) + str(i)
        df_i = pd.read_csv(i_dir,sep="\t")

        # Create unique ID for each individual.
        df_i['ID'] = df_i['#FID'].astype(str) + "_" + df_i['IID'].astype(str)
        
        # Extract only ID and score sum column
        df_i = df_i[['ID','SCORE1_SUM']]
        
        # Rename SCORE1_SUM column so that it is unique to each chr.
        df_i.columns = ['ID',i]

        # Merge all chr scores to df.
        if len(merged_df) == 0:
            merged_df = df_i
        else:
            merged_df = pd.merge(df_i,merged_df,on="ID",how="outer")

    # Sum the scores across chr for each individual.
    merged_df['PRS'] = merged_df.loc[:,score_lst].sum(axis=1)

    # Extract individual ID and final PRS calculation.
    merged_df = merged_df[['ID','PRS']]

    # Write final PRS calculations to file.
    merged_df.to_csv(out_file_name,sep="\t",index=False,header=True)
    
    print(merged_df)

    return True 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", help = "Please input the path to the directory with the final .sscore files for final PRS generation.")
    parser.add_argument("-o", help = "Please input the path to the directory you would like to put your final PRS calculation in.")
    args = parser.parse_args()
    
    scoreLst = score_dir(args.d)

    computePRS = compute_PRS(args.d,scoreLst,args.o)
    
if __name__ == "__main__":
    main()
