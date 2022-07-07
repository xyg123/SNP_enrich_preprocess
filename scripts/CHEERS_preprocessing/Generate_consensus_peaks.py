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

def main():

    start = time.time()
    print("Generating consensus peak file...")
    args = parse_args()
    #for 
    Input_peaks=pd.read_csv(args.Peaks, sep="\t")
    # Current peaks =
    for i in range(0, len(Input_peaks.index)):
        
        urllib.request.urlretrieve(Input_peaks.iloc[i]["URL"], "../tmp/tmp_peak.bed.gz")
        Current_peaks=pybedtools.BedTool("../tmp/tmp_peak.bed.gz")
        Current_peaks_pd=pd.read_table(Current_peaks.fn)
        Current_peaks_pd.to_csv("../tmp/Concat_peaks.bed", mode='a', header=False, sep="\t", index=False)

    # Formatted peaks =

    print(str(len(Input_peaks.index))+" peak files read in "+str(time.time()-start)+" seconds...")
    
    Concat_peaks=pybedtools.BedTool("../tmp/Concat_peaks.bed")
    Concat_peaks_sorted=Concat_peaks.sort()
    Concat_peaks_merged=Concat_peaks_sorted.merge(d=0)

    Concat_peaks_merged.saveas(args.outdir+args.prefix+"_Consensus_peaks.bed")

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', metavar="<str>", help=("Output prefix"), type=str, required=True)
    parser.add_argument('--Peaks', metavar="<str>", help=("Input peak URLs"), type=str, required=True)
    parser.add_argument('--outdir', metavar="<str>", help=("Output directory"), type=str, required=True)
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()

# python Generate_consensus_peaks.py --prefix BLUEPRINT --Peaks ~/BLUEPRINT_peak_URLs.tsv --outdir ~/BLUEPRINT_peaks/
#  
