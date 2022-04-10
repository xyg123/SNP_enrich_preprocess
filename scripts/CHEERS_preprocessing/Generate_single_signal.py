import os
import sys
import numpy as np
import pandas as pd
import argparse
import pybedtools
import re
import pyBigWig as pbw
import time
import urllib.request




def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--Sample', metavar="<str>", help=("Input Sample ID"), type=str, required=True)
    parser.add_argument('--Peaks', metavar="<str>", help=("Input Consensus peaks"), type=str, required=True)
    parser.add_argument('--outdir', metavar="<str>", help=("Output"), type=str, required=True)
    
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    start = time.time()
    print("Generating signal for "+args.Sample+" consensus peak file...")

    imputed_urls_path="../../configs/Imputed_urls.txt"
    observed_urls_path="../../configs/Observed_urls.txt"

    Imputed_urls=pd.read_csv(imputed_urls_path, header=None, names=['filename'])
    Observed_urls=pd.read_csv(observed_urls_path, header=None, names=['filename'])

    sample_ID=args.Sample
    Master_merged_pd=pd.read_table(args.Peaks, 
                            header=None, names=['chr', 'start', 'end'])
                            
    BW_URL=Get_BW_URLs(sample_ID, "H3K27ac", Imputed_urls, Observed_urls)
    print("Bigwig URL obtained ...")
    Current_BW=pbw.open(BW_URL)
    print("Generating signals at peaks now ...")
    Signals=Master_merged_pd.apply(lambda x: sum(Current_BW.values(x['chr'], x['start'], x['end'])), axis=1)
    Master_merged_pd['signal']=Signals

    Master_merged_pd.to_csv(args.outdir+"/"+sample_ID+"_ReadsInPeaks.txt",
                        sep="\t", index=False, header=['chrom', 'start', 'end', sample_ID])

    end = time.time()
    
    print(end - start)

def Get_BW_URLs(Sample_ID, Chromatin_mark, imputed_urls, observed_urls):
    Imputed_search=imputed_urls.loc[imputed_urls['filename'].str.contains(Sample_ID+"_"+Chromatin_mark)].reset_index()
    if Imputed_search.empty == False:
        Matched_URLs="https://epigenome.wustl.edu/epimap/data/imputed/"+Imputed_search.iloc[0]['filename']
    Observed_search=observed_urls.loc[observed_urls['filename'].str.contains(Chromatin_mark+"_"+Sample_ID)].reset_index()
    if Observed_search.empty == False:
        Matched_URLs="https://epigenome.wustl.edu/epimap/data/observed/"+Observed_search.iloc[0]['filename']
    return(Matched_URLs)

if __name__ == '__main__':
    
    main()

# python Generate_single_signal.py --Sample "BSS01668" --Peaks "../../tmp/Master_enhancers_d250.sorted.merged.bed" --outdir "../../tmp/H3K27ac" 