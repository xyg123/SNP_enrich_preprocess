import os
import sys
import numpy as np
import pandas as pd

import pybedtools
import re
import pyBigWig as pbw
import time
import urllib.request
import matplotlib.pyplot as plt

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

def strip_to_sample(input_url):
    return(input_url.split(sep="_")[0])

Promoter_urls_samples=Promoter_urls.apply(lambda x: strip_to_sample(x['filename']), axis=1)
Enhancer_urls_samples=Enhancer_urls.apply(lambda x: strip_to_sample(x['filename']), axis=1)

epimap_meta=epimap_meta.loc[epimap_meta['id'].isin(Enhancer_urls_samples) & epimap_meta['id'].isin(Promoter_urls_samples)]

Test_samples=["BSS00121", "BSS00122", "BSS00123", "BSS00206", "BSS00207", "BSS01712", "BSS01714", "BSS00196", "BSS00197", "BSS00188", "BSS00189", "BSS00233", "BSS00234", "BSS00227", "BSS00228", "BSS01113", "BSS01120", "BSS00084", "BSS00330", "BSS00079", "BSS00080", "BSS01665", "BSS01666", "BSS01169", "BSS00511", "BSS01068", "BSS01071"]

epimap_meta=epimap_meta.loc[epimap_meta['id'].isin(Test_samples)]
#Cell_types_to_run=['Digestive', 'Blood & T-cell', 'Bone', 'Brain', 'Adipose', 'Mesench', 
#    'Endocrine', 'Stromal', 'Epithelial', 'Heart', 'HSC & B-cell', 'ES-deriv', 'Pancreas',
#    'Endothelial'];

Cell_types_to_run=epimap_meta['GROUP'].unique()
#Cell_types_to_run=np.delete(Cell_types_to_run, range(0, 7))



def main():
    start = time.time()
    print("Generating consensus peak file...")
    for cell in Cell_types_to_run:
        cell_IDs=Get_Sample_IDs(cell)
        for sample_ID in cell_IDs:
            Reg_URLs=Get_regulatory_URLs(sample_ID)
            urllib.request.urlretrieve(Reg_URLs[0], "../tmp/temp_enhancer.bed.gz")
            urllib.request.urlretrieve(Reg_URLs[1], "../tmp/temp_promoter.bed.gz")
            
            Current_enhancers=pybedtools.BedTool("../tmp/temp_enhancer.bed.gz")
            Current_enhancers_pd=pd.read_table(Current_enhancers.fn)
            Current_enhancers_pd.to_csv("../output/Master_enhancers.bed",
                                        mode='a', header=False, sep="\t", index=False)
            
            Current_promoters=pybedtools.BedTool("../tmp/temp_promoter.bed.gz")
            Current_promoters_pd=pd.read_table(Current_promoters.fn)
            Current_promoters_pd.to_csv("../output/Master_enhancers.bed",
                                        mode='a', header=False, sep="\t", index=False)
            
    # Merge the Promoter+enhancer files

    Master_enhancers=pybedtools.BedTool("../output/Master_enhancers.bed")
    Master_enhancers_sorted=Master_enhancers.sort()
    Master_enhancers_merged=Master_enhancers_sorted.merge(d=0)

    Master_enhancers_merged.saveas("../output/Master_enhancers.sorted.merged.bed")

    Master_merged_pd=pd.read_table("../output/Master_enhancers.sorted.merged.bed", 
                              header=None, names=['chr', 'start', 'end'])
                              
    print("Generating signal files...")
# Loop through bigwig URLs and get signal counts in merged peaks

    for cell in Cell_types_to_run:
        cell_IDs=Get_Sample_IDs(cell)
        for sample_ID in cell_IDs:
            BW_URL=Get_BW_URLs(sample_ID, "H3K27ac")
            Current_BW=pbw.open(BW_URL)
            Signals=Master_merged_pd.apply(lambda x: sum(Current_BW.values(x['chr'], x['start'], x['end'])), axis=1)
            Master_merged_pd['signal']=Signals

            Master_merged_pd.to_csv("../output/H3K27ac_ReadsInPeaks/"+sample_ID+"_ReadsInPeaks.txt",
                                sep="\t", index=False, header=['chrom', 'start', 'end', sample_ID])

    end = time.time()
    print(end - start)


def Get_Sample_IDs(Cell_type):
    Matched_celltype=epimap_meta.loc[epimap_meta['GROUP'].str.contains(Cell_type)].reset_index()
    return(Matched_celltype['id'])
    

def Get_regulatory_URLs(Sample_ID):
    enhancer_page="https://personal.broadinstitute.org/cboix/epimap/mark_matrices/enhancers_bysample/"
    promoter_page="https://personal.broadinstitute.org/cboix/epimap/mark_matrices/promoters_bysample/"
    
    enhancer_file=Enhancer_urls.loc[Enhancer_urls['filename'].str.contains(Sample_ID)].reset_index()
    promoter_file=Promoter_urls.loc[Promoter_urls['filename'].str.contains(Sample_ID)].reset_index()
    
    return([enhancer_page+enhancer_file.iloc[0]['filename'], 
            promoter_page+promoter_file.iloc[0]['filename']])

def Get_BW_URLs(Sample_ID, Chromatin_mark):
    Imputed_search=imputed_urls.loc[imputed_urls['filename'].str.contains(Sample_ID+"_"+Chromatin_mark)].reset_index()
    if Imputed_search.empty == False:
        Matched_URLs="https://epigenome.wustl.edu/epimap/data/imputed/"+Imputed_search.iloc[0]['filename']
    Observed_search=observed_urls.loc[observed_urls['filename'].str.contains(Chromatin_mark+"_"+Sample_ID)].reset_index()
    if Observed_search.empty == False:
        Matched_URLs="https://epigenome.wustl.edu/epimap/data/observed/"+Observed_search.iloc[0]['filename']
    return(Matched_URLs)


if __name__ == '__main__':

    main()


