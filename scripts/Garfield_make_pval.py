import pandas as pd
import os
import sys
import numpy as np
import argparse
from subprocess import Popen, PIPE
import gzip
from pybedtools import BedTool
from collections import OrderedDict
import time
from pyliftover import LiftOver

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

def Garfield_make_pval_with_LO(input_sumstats, trait, out_path, lo_from, lo_to):
    start_time=time.time()
    lo=LiftOver(lo_from, lo_to)
    Sumstats=pd.read_parquet(input_sumstats, engine="pyarrow")
    Sumstats['chrom']="chr"+Sumstats['chrom']
    Sumstats['lift_chr']=Sumstats.apply(lambda x: Liftover_chr(input_chr=x['chrom'], input_pos=x['pos'], input_LO=lo), axis=1)
    Sumstats['lift_pos']=Sumstats.apply(lambda x: Liftover_pos(input_chr=x['chrom'], input_pos=x['pos'], input_LO=lo), axis=1)
    Sumstats=Sumstats.loc[-Sumstats['lift_chr'].str.contains("NA")]
    cols = OrderedDict([
        ('lift_chr', 'chrom'),
        ('lift_pos', 'pos'),
        ('pval', 'pval')
        ])
    Sumstats=Sumstats.loc[:, list(cols.keys())].rename(columns=cols) 
    print("finished liftover of "+lo_from+" to "+lo_to+" in "+str(time.time()-start_time)+" seconds!")

    for i in range(1, 23):
        start_time=time.time()
        Sumstats_bychr=Sumstats.loc[Sumstats['chrom'] == str(i)][["pos", "pval"]]
        Sumstats_bychr.to_csv(out_path+trait+"/chr"+str(i), sep="\t", index=False, header=None)
        print("finished chr"+str(i)+" in "+str(time.time()-start_time)+" seconds!")
    start_time=time.time()
    Sumstats_bychr=Sumstats.loc[Sumstats['chrom'] == "X"][["pos", "pval"]]
    Sumstats_bychr.to_csv(out_path+trait+"/chrX", sep="\t", index=False, header=None)
    print("finished chrX"+" in "+str(time.time()-start_time)+" seconds!")

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('--Input_sumstats', metavar="<str>", type=str, required=True)
    parser.add_argument('--Trait', metavar="<str>", type=str, required=True)
    parser.add_argument('--Output_dir', metavar="<str>", type=str, required=True)
    parser.add_argument('--Lift_from', metavar="<str>", type=str, required=True)
    parser.add_argument('--Lift_to', metavar="<str>", type=str, required=True)
    
    args = parser.parse_args()
    return args

def main():
    args=parse_args()
    Garfield_make_pval_with_LO(input_sumstats=args.Input_sumstats, out_path=args.Output_dir, trait=args.Trait, lo_from=args.Lift_from, lo_to=args.Lift_to)

if __name__ == '__main__':

    main()
