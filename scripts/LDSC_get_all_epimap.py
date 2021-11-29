import pandas as pd
import os
import sys

import numpy as np

import argparse
import re

import gzip
from pybedtools import BedTool

from collections import OrderedDict
import time

import urllib.request

import json
import requests

def main():

    args = parse_args()

    imputed_urls_path="../configs/Imputed_urls.txt"
    observed_urls_path="../configs/Observed_urls.txt"

    Enhancer_urls=pd.read_csv("../configs/enhancers_urls.txt",
                            header=None, names=['filename'])
    Promoter_urls=pd.read_csv("../configs/promoters_urls.txt",
                            header=None, names=['filename'])

    imputed_urls=pd.read_csv(imputed_urls_path, header=None, names=['filename'])
    observed_urls=pd.read_csv(observed_urls_path, header=None, names=['filename'])

    epimap_meta_path="../configs/main_metadata_table.tsv"

    epimap_meta=pd.read_csv(epimap_meta_path, sep="\t")

    Promoter_urls_samples=Promoter_urls.apply(lambda x: strip_to_sample(x['filename']), axis=1)
    Enhancer_urls_samples=Enhancer_urls.apply(lambda x: strip_to_sample(x['filename']), axis=1)

    epimap_meta=epimap_meta.loc[epimap_meta['id'].isin(Enhancer_urls_samples) & epimap_meta['id'].isin(Promoter_urls_samples)]

    Cell_types_to_run=epimap_meta['GROUP'].unique()

    for cell in Cell_types_to_run:
        cell_IDs=Get_Sample_IDs(cell, epimap_meta)
        for sample_ID in cell_IDs:
            Reg_URLs=Get_regulatory_URLs(sample_ID, Enhancer_urls, Promoter_urls)
            urllib.request.urlretrieve(Reg_URLs[0], args.outdir+sample_ID+"_enhancer.bed.gz")
            urllib.request.urlretrieve(Reg_URLs[1], args.outdir+sample_ID+"_promoter.bed.gz")


def strip_to_sample(input_url):
    return(input_url.split(sep="_")[0])



def Get_Sample_IDs(Cell_type, Meta_file):
    Matched_celltype=Meta_file.loc[Meta_file['GROUP'].str.contains(Cell_type)].reset_index()
    return(Matched_celltype['id'])
    

def Get_regulatory_URLs(Sample_ID, E_URLs, P_URLs):
    enhancer_page="https://personal.broadinstitute.org/cboix/epimap/mark_matrices/enhancers_bysample/"
    promoter_page="https://personal.broadinstitute.org/cboix/epimap/mark_matrices/promoters_bysample/"
    
    enhancer_file=E_URLs.loc[E_URLs['filename'].str.contains(Sample_ID)].reset_index(drop=True)
    promoter_file=P_URLs.loc[P_URLs['filename'].str.contains(Sample_ID)].reset_index(drop=True)
    
    return([enhancer_page+enhancer_file.iloc[0]['filename'], 
            promoter_page+promoter_file.iloc[0]['filename']])

# Get their associated Promoter+Enhancer files:



def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', metavar="<str>", help=("Output"), type=str, required=True)
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()