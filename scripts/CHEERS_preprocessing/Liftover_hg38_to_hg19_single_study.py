import os
import sys
import json
import argparse
import gzip
import pandas as pd
import pyarrow
from pyliftover import LiftOver


# Gathers all credible SNPs from the GCS finemapping file, liftover from hg38 to hg19,
#  then outputs a run command for CHEERS : 
#  i.e. python3 CHEERS_computeEnrichment.py 
#           Height (Use Study ID)
#           Results/ (Use Enrichment_outdir)
#           (Use input_signal) ../Subsampled_epimap_H3K27ac/Subsampled_Epimap_H3K27ac_counts_normToMax_quantileNorm_euclideanNorm.txt 
#           (Use outdir+Study_ID+"_ hg19.txt") GWAS_inputs/NEALE2_50_raw_Height_hg19.txt

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--Study_ID', metavar="<str>", help=("Input study ID"), type=str, required=True)
    parser.add_argument('--input_credset', metavar="<str>", help=("Input Credible SNP set file"), type=str, required=True)
    parser.add_argument('--Enrichment_outdir', metavar="<str>", help=("Output"), type=str, required=True)
    parser.add_argument('--input_peak', metavar="<str>", help=("Input signal at peaks file"), type=str, required=True)
    parser.add_argument('--outdir', metavar="<str>", help=("Output for hg19 SNPs"), type=str, required=True)
    
    args = parser.parse_args()
    return args

def Call_liftover(input_chr, input_pos, input_LO):
    return(input_LO.convert_coordinate(input_chr, input_pos))

def Liftover_chr(input_chr, input_pos, input_LO):
    lifted_result=input_LO.convert_coordinate(input_chr, input_pos)
    if(len(lifted_result) == 0):
        return("NA")
    else:
        return(lifted_result[0][0])

def Liftover_pos(input_chr, input_pos, input_LO):
    lifted_result=input_LO.convert_coordinate(input_chr, input_pos)
    if(len(lifted_result) == 0):
        return("NA")
    else:
        return(lifted_result[0][1])

def main():
    args=parse_args()
    lo=LiftOver("hg38", "hg19")

    All_credsets=pd.read_parquet(args.input_credset, engine='pyarrow')

    Study_SNPs=All_credsets.loc[All_credsets["study_id"].str.contains(args.Study_ID)].reset_index(drop=True)
    Study_SNPs['tag_chrom']="chr"+Study_SNPs['tag_chrom']
    Study_SNPs['hg19_chr']=Study_SNPs.apply(lambda x: Liftover_chr(input_chr=x['tag_chrom'], input_pos=x['tag_pos'], input_LO=lo), axis=1)
    Study_SNPs['hg19_pos']=Study_SNPs.apply(lambda x: Liftover_pos(input_chr=x['tag_chrom'], input_pos=x['tag_pos'], input_LO=lo), axis=1)
    Study_SNPs=Study_SNPs.loc[-Study_SNPs['hg19_chr'].str.contains("NA")]

    #Study_SNPs['hg19_names']=Study_SNPs.apply(lambda x: x['tag_chrom']+"_"+str(x['tag_pos']), axis=1)
    #Study_SNPs['tag_chrom']="chr"+Study_SNPs['tag_chrom']
    #Study_SNPs=Study_SNPs[['hg19_names', 'tag_chrom', 'tag_pos']]
    
    Study_SNPs['hg19_names']=Study_SNPs.apply(lambda x: x['hg19_chr']+"_"+str(x['hg19_pos']), axis=1)
    Study_SNPs=Study_SNPs[['hg19_names', 'hg19_chr', 'hg19_pos']]
    Study_SNPs.to_csv(args.outdir+args.Study_ID+"_hg19.txt", sep="\t", index=False, header=False)

    cmd=['python3', "CHEERS_computeEnrichment.py", args.Study_ID, args.Enrichment_outdir, args.input_peak, args.outdir+args.Study_ID+"_hg19.txt"]
    cmd_str=" ".join([str(arg) for arg in cmd])
    print(cmd_str)
    out_todo="~/compute_CHEERS_enrichments.txt"
    todo_h=open(os.path.expanduser(out_todo), "ab")
    todo_h.write((cmd_str + '\n').encode())

if __name__ == '__main__':
    
    main()

# python Liftover_hg38_to_hg19_single_study.py 
#   --Study_ID GCSTxxxx 
#   --input_credset ~/Credible_SNP_sets/finemapping_220401.parquet/
#   --Enrichment_outdir ~/CHEERS/Results/
#   --input_peak ~/Subsampled_epimap_H3K27ac/Subsampled_Epimap_H3K27ac_counts_normToMax_quantileNorm_euclideanNorm.txt
#   --outdir ~/Credible_SNP_sets/Formatted_hg19/
