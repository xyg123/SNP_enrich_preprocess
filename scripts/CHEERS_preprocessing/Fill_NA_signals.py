import os
import sys
import numpy as np
import pandas as pd
import argparse

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--File', metavar="<str>", help=("Input Sample ID"), type=str, required=True)

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    Signals=pd.read_table(args.File, sep="\t")
    Sample=args.File.split("_")[0]

    Signals[Sample]=Signals[Sample].fillna(0)

    Signals.to_csv(args.File, sep="\t", index=False)


if __name__ == '__main__':
    
    main()
