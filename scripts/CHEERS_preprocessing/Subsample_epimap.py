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

imputed_urls_path="../../configs/Imputed_urls.txt"
observed_urls_path="../../configs/Observed_urls.txt"

Enhancer_urls=pd.read_csv("../../configs/enhancers_urls.txt",
                        header=None, names=['filename'])
Promoter_urls=pd.read_csv("../../configs/promoters_urls.txt",
                        header=None, names=['filename'])

imputed_urls=pd.read_csv(imputed_urls_path, header=None, names=['filename'])
observed_urls=pd.read_csv(observed_urls_path, header=None, names=['filename'])

epimap_meta_path="../../configs/main_metadata_table.tsv"

epimap_meta=pd.read_csv(epimap_meta_path, sep="\t")

def strip_to_sample(input_url):
    return(input_url.split(sep="_")[0])

Promoter_urls_samples=Promoter_urls.apply(lambda x: strip_to_sample(x['filename']), axis=1)
Enhancer_urls_samples=Enhancer_urls.apply(lambda x: strip_to_sample(x['filename']), axis=1)

epimap_meta=epimap_meta.loc[epimap_meta['id'].isin(Enhancer_urls_samples) & epimap_meta['id'].isin(Promoter_urls_samples)]

#Test_samples=["BSS00121", "BSS00122", "BSS00123", "BSS00206", "BSS00207", "BSS01712", "BSS01714", "BSS00196", "BSS00197", "BSS00188", "BSS00189", "BSS00233", "BSS00234", "BSS00227", "BSS00228", "BSS01113", "BSS01120", "BSS00084", "BSS00330", "BSS00079", "BSS00080", "BSS01665", "BSS01666", "BSS01169", "BSS00511", "BSS01068", "BSS01071"]
Test_samples=['BSS01668', 'BSS01665', 'BSS00038', 'BSS01666', 'BSS01669', 'BSS01696', 'BSS01690', 'BSS00183', 'BSS01693', 'BSS01420', 'BSS01154', 'BSS00330', 'BSS00084', 'BSS01397', 'BSS00705', 'BSS00078', 'BSS00077', 'BSS00131', 'BSS01272', 'BSS00127', 'BSS00025', 'BSS00748', 'BSS00026', 'BSS00159', 'BSS01105', 'BSS01548', 'BSS01594', 'BSS01657', 'BSS01282', 'BSS01602', 'BSS00273', 'BSS01264', 'BSS00287', 'BSS01612', 'BSS01366', 'BSS00716', 'BSS00277', 'BSS00478', 'BSS00315', 'BSS00483', 'BSS00057', 'BSS00050', 'BSS00284', 'BSS01399', 'BSS00048', 'BSS00258', 'BSS01206', 'BSS00143', 'BSS00260', 'BSS00387', 'BSS01505', 'BSS00355', 'BSS00075', 'BSS01068', 'BSS01103', 'BSS01503', 'BSS00329', 'BSS01504', 'BSS00328', 'BSS01502', 'BSS00546', 'BSS00238', 'BSS00098', 'BSS00550', 'BSS00237', 'BSS00495', 'BSS01127', 'BSS00087', 'BSS01508', 'BSS00520', 'BSS01078', 'BSS01136', 'BSS01513', 'BSS01482', 'BSS01533', 'BSS00554', 'BSS01158', 'BSS01170', 'BSS00553', 'BSS00511', 'BSS01871', 'BSS01527', 'BSS01147', 'BSS01187', 'BSS01188', 'BSS00428', 'BSS00452', 'BSS00404', 'BSS00439', 'BSS00472', 'BSS01315', 'BSS01301', 'BSS01846', 'BSS01324', 'BSS01298', 'BSS01574', 'BSS01573', 'BSS01338', 'BSS01155', 'BSS01576', 'BSS00146', 'BSS00145', 'BSS00304', 'BSS00368', 'BSS01156', 'BSS01621', 'BSS01840', 'BSS01617', 'BSS01842', 'BSS01613', 'BSS00123', 'BSS00122', 'BSS00124', 'BSS00121', 'BSS00758', 'BSS01855', 'BSS01856', 'BSS01443', 'BSS01859', 'BSS01867', 'BSS01456', 'BSS01886', 'BSS01884', 'BSS01887', 'BSS01459', 'BSS01286', 'BSS01660', 'BSS01606', 'BSS01287', 'BSS01659', 'BSS01625', 'BSS01628', 'BSS01634', 'BSS01631', 'BSS01630', 'BSS00339', 'BSS00332', 'BSS01661', 'BSS01583', 'BSS00063', 'BSS01823', 'BSS01827', 'BSS01824', 'BSS01825', 'BSS01826', 'BSS00244', 'BSS00734', 'BSS01108', 'BSS01107', 'BSS00738']

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


def main():
    start = time.time()
    print("Generating consensus peak file...")
    for sample_ID in Test_samples:
        Reg_URLs=Get_regulatory_URLs(sample_ID)
        urllib.request.urlretrieve(Reg_URLs[0], "../../tmp/temp_enhancer.bed.gz")
        urllib.request.urlretrieve(Reg_URLs[1], "../../tmp/temp_promoter.bed.gz")
        
        Current_enhancers=pybedtools.BedTool("../../tmp/temp_enhancer.bed.gz")
        Current_enhancers_pd=pd.read_table(Current_enhancers.fn)
        Current_enhancers_pd.to_csv("../../tmp/Master_enhancers.bed",
                                    mode='a', header=False, sep="\t", index=False)
        
        Current_promoters=pybedtools.BedTool("../../tmp/temp_promoter.bed.gz")
        Current_promoters_pd=pd.read_table(Current_promoters.fn)
        Current_promoters_pd.to_csv("../../tmp/Master_enhancers.bed",
                                    mode='a', header=False, sep="\t", index=False)
            
    # Merge the Promoter+enhancer files

    Master_enhancers=pybedtools.BedTool("../../tmp/Master_enhancers.bed")
    Master_enhancers_sorted=Master_enhancers.sort()
    Master_enhancers_merged=Master_enhancers_sorted.merge(d=0)

    Master_enhancers_merged.saveas("../../tmp/Master_enhancers.sorted.merged.bed")

    Master_merged_pd=pd.read_table("../../tmp/Master_enhancers.sorted.merged.bed", 
                              header=None, names=['chr', 'start', 'end'])
                              
    print("Generating signal files...")
# Loop through bigwig URLs and get signal counts in merged peaks

    for sample_ID in Test_samples:
        BW_URL=Get_BW_URLs(sample_ID, "H3K27ac")
        Current_BW=pbw.open(BW_URL)
        Signals=Master_merged_pd.apply(lambda x: sum(Current_BW.values(x['chr'], x['start'], x['end'])), axis=1)
        Master_merged_pd['signal']=Signals

        Master_merged_pd.to_csv("../../tmp/H3K27ac_ReadsInPeaks/"+sample_ID+"_ReadsInPeaks.txt",
                            sep="\t", index=False, header=['chrom', 'start', 'end', sample_ID])

    end = time.time()
    print(end - start)

if __name__ == '__main__':
    
    main()

