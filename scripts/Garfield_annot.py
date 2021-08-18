import pandas as pd
import os
import sys
#import fastparquet
import numpy as np
#import tabix

from subprocess import Popen, PIPE

import gzip
from pybedtools import BedTool

from collections import OrderedDict
import time

def Garfield_annot_UK10K(input_SNP_path, out_path, bed_for_annot_dir):
    for i in range(1, 23):
        start_time=time.time()

        Input_SNPs=pd.read_csv(input_SNP_path+"chr"+str(i), names=["BP", "MAF", "TSS"], sep=" ")
        Input_SNPs['Chrom']="chr"+str(i)

        To_overlap=Input_SNPs[['Chrom', 'BP']]
        To_overlap['start']=To_overlap['BP'].astype(int)-1
        To_overlap=To_overlap[['Chrom', 'start', 'BP']]

        Annot_out=To_overlap['BP']

        SNP_bed=BedTool.from_dataframe(To_overlap)
        for bed_for_annot in os.listdir(bed_for_annot_dir):
            Annot_bed=BedTool(bed_for_annot_dir+bed_for_annot)

            SNP_intersect=SNP_bed.intersect(Annot_bed)
            bp=[x.start + 1 for x in SNP_intersect]
            df_int=pd.DataFrame({'BP': bp, bed_for_annot: int(1)})
            Annot_out=pd.merge(Annot_out, df_int, how="left", on='BP')
            Annot_out.fillna(0, inplace=True)

        Annot_out=Annot_out.astype(int)
        Annot_out['annotations']=Annot_out.iloc[:,1:].apply(lambda row: "".join(row.values.astype(str)), axis=1)
        Annot_out=Annot_out[['BP', 'annotations']]
        Annot_out.to_csv(out_path+"chr"+str(i), sep=" ", header=None, index=False)
        print("finished chr"+str(i)+" in "+str(time.time()-start_time)+" seconds!")
        #return(Annot_out)

#Garfield_annot_UK10K("/Users/xg1/Downloads/garfield-data/maftssd/", "/Users/xg1/Downloads/garfield-data/Comparison_annotations/", "/Users/xg1/Downloads/ldsc/Comparison_enhancers/")

def Garfield_annot_UK10K_line(input_SNP_path, out_path, bed_for_annot_dir):
    for i in range(1, 23):
        start_time=time.time()
        with open(input_SNP_path+"chr"+str(i), 'r') as f:
            for line in f:
                line=line.rstrip().split(" ")
                bp=line[0];
                input_chr="chr"+str(i)
                input_start=int(bp)-1
                input_SNP=BedTool(str(input_chr)+" "+str(input_start)+" "+str(bp), from_string=True);

                annot_out=[]

                for bed_for_annot in os.listdir(bed_for_annot_dir):
                    Annot_bed=BedTool(bed_for_annot_dir+bed_for_annot)
                    SNP_intersect=input_SNP.intersect(Annot_bed)
                    annot_out.append(str(len(SNP_intersect)))
                
                with open(out_path+"chr"+str(i), 'a') as out_file:
                    out_file.write(str(bp)+" "+"".join(annot_out)+"\n")

        # Input_SNPs=pd.read_csv(input_SNP_path+"chr"+str(i), names=["BP", "MAF", "TSS"], sep=" ")
        # Input_SNPs['Chrom']="chr"+str(i)

        # To_overlap=Input_SNPs[['Chrom', 'BP']]
        # To_overlap['start']=To_overlap['BP'].astype(int)-1
        # To_overlap=To_overlap[['Chrom', 'start', 'BP']]

        # Annot_out=To_overlap['BP']

        # SNP_bed=BedTool.from_dataframe(To_overlap)
        # for bed_for_annot in os.listdir(bed_for_annot_dir):
        #     Annot_bed=BedTool(bed_for_annot_dir+bed_for_annot)

        #     SNP_intersect=SNP_bed.intersect(Annot_bed)
        #     bp=[x.start + 1 for x in SNP_intersect]
        #     df_int=pd.DataFrame({'BP': bp, bed_for_annot: int(1)})
        #     Annot_out=pd.merge(Annot_out, df_int, how="left", on='BP')
        #     Annot_out.fillna(0, inplace=True)

        # Annot_out=Annot_out.astype(int)
        # Annot_out['annotations']=Annot_out.iloc[:,1:].apply(lambda row: "".join(row.values.astype(str)), axis=1)
        # Annot_out=Annot_out[['BP', 'annotations']]
        # Annot_out.to_csv(out_path+"chr"+str(i), sep=" ", header=None, index=False)
        print("finished chr"+str(i)+" in "+str(time.time()-start_time)+" seconds!")
        #return(annot_out)
        #return(Annot_out)


abc=Garfield_annot_UK10K_line("/Users/xg1/Downloads/garfield-data/maftssd/", 
    "/Users/xg1/Downloads/garfield-data/Comparison_annotations/",
    "/Users/xg1/Downloads/ldsc/Comparison_enhancers/")