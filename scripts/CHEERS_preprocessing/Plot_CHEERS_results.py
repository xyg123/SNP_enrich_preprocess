import os
import sys
import numpy as np
import pandas as pd
#from scipy.stats import norm

import pybedtools
import re
#import matplotlib as plt
import pyBigWig as pbw
import time
import urllib.request
import matplotlib.pyplot as plt

epimap_meta_path="../configs/main_metadata_table.tsv"
epimap_meta=pd.read_csv(epimap_meta_path, sep="\t")

def add_GROUP(CHEERS_res):
    Sample_group=[]
    Sample_info=[]
    for i in CHEERS_res['Sample']:
        tmp=epimap_meta.loc[epimap_meta['id'].str.contains(i)]
        #print(tmp.iloc[0]['GROUP'])
        Sample_group.append(tmp.iloc[0]['GROUP'])
        Sample_info.append(tmp.iloc[0]['infoline'])
    CHEERS_res['GROUP']=Sample_group
    CHEERS_res['infoline']=Sample_info
    return(CHEERS_res)

Corrected_IBD_P=pd.read_csv("/home/xg1/CHEERS/Results/IBD_GCST_DNase_disease_enrichment_pValues.txt",
                      sep="\t", header=None, names=["Sample", "P"])
Corrected_IBD_P=add_GROUP(Corrected_IBD_P)
print(Corrected_IBD_P)

