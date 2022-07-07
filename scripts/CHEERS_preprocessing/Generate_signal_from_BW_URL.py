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
    parser.add_argument('--BW_URL', metavar="<str>", help=("Input bigWig URL"), type=str, required=True)
    parser.add_argument('--outdir', metavar="<str>", help=("Output"), type=str, required=True)
    
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    start = time.time()
    print("Generating signal for "+args.Sample+" consensus peak file...")

    Master_merged_pd=pd.read_table(args.Peaks, 
                            header=None, names=['chr', 'start', 'end'])

    Current_BW=pbw.open(args.BW_URL)
    print("Generating signals at peaks now ...")
    Signals=Master_merged_pd.apply(lambda x: sum(Current_BW.values(x['chr'], x['start'], x['end'])), axis=1)
    Master_merged_pd['signal']=Signals

    Master_merged_pd.to_csv(args.outdir+"/"+args.Sample+"_ReadsInPeaks.txt",
                        sep="\t", index=False, header=['chrom', 'start', 'end', args.Sample])

    end = time.time()
    
    print(end - start)

if __name__ == '__main__':
    
    main()

# python Generate_single_signal.py --Sample "BSS01668" --Peaks "../../tmp/Master_enhancers_d250.sorted.merged.bed" --outdir "../../tmp/H3K27ac" 