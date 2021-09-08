import pandas as pd
import os
import sys
#import fastparquet
import numpy as np
#import tabix
import argparse
import re
from subprocess import Popen, PIPE

import gzip
from pybedtools import BedTool

from collections import OrderedDict
import time

epimap_meta_path="../configs/main_metadata_table.tsv"

def add_GROUP(Input_sample, meta_file):
    epimap_meta=pd.read_csv(meta_file, sep="\t")
    tmp=epimap_meta.loc[epimap_meta['id'].str.contains(Input_sample)]
    if(len(tmp) > 0):
        return(tmp['GROUP'])
    else:
        return("NA")

def add_infoline(Input_sample, meta_file):
    epimap_meta=pd.read_csv(meta_file, sep="\t")
    tmp=epimap_meta.loc[epimap_meta['id'].str.contains(Input_sample)]
    if(len(tmp) > 0):
        return(tmp['infoline'])
    else:
        return("NA")
  

def test_epimap(Input_sample, meta_file):
    epimap_meta=pd.read_csv(meta_file, sep="\t")
    tmp=epimap_meta.loc[epimap_meta['id'].str.contains(Input_sample)]
    return(tmp)

def Garfield_annot_UK10K(input_SNP_path, out_path, bed_for_annot_dir):
    for i in range(1, 23):
        start_time=time.time()
        if not os.path.exists(out_path+"chr"+str(i)):
            print("{} does not exist, creating".format(out_path+"chr"+str(i)))
            os.mkdir(out_path+"chr"+str(i)+"/")

        Input_SNPs=pd.read_csv(input_SNP_path+"chr"+str(i), names=["BP", "MAF", "TSS"], sep=" ")
        Input_SNPs['Chrom']="chr"+str(i)
        Input_SNPs=Input_SNPs.drop(["MAF", "TSS"], axis=1)

        Input_SNPs['start']=Input_SNPs['BP'].astype(int)-1
        Input_SNPs=Input_SNPs[['Chrom', 'start', 'BP']]

        #Annot_out=Input_SNPs['BP']

        SNP_bed=BedTool.from_dataframe(Input_SNPs)
        for bed_for_annot in os.listdir(bed_for_annot_dir):
            Annot_out=Input_SNPs['BP']
            start_time_chr=time.time()

            Annot_bed=BedTool(bed_for_annot_dir+bed_for_annot)
            Annot_bed=Annot_bed.merge()
            Sample_name=bed_for_annot.split(".")[0]

            SNP_intersect=SNP_bed.intersect(Annot_bed)
            bp=[x.start + 1 for x in SNP_intersect]
            df_int=pd.DataFrame({'BP': bp, bed_for_annot: int(1)})
            Annot_out=pd.merge(Annot_out, df_int, how="left", on='BP')
            Annot_out.fillna(0, inplace=True)

            Annot_out=Annot_out.astype(int)
            
            Annot_out.to_csv(out_path+"chr"+str(i)+"/"+Sample_name, sep=" ", index=False)
            print("finished chr"+str(i)+" for "+Sample_name+" in "+str(time.time()-start_time_chr)+" seconds!")
            #print("finished chr"+str(i)+" in "+str(time.time()-start_time_chr)+" seconds!")

    Input_SNPs=pd.read_csv(input_SNP_path+"chrX", names=["BP", "MAF", "TSS"], sep=" ")
    Input_SNPs['Chrom']="chrX"
    Input_SNPs=Input_SNPs.drop(["MAF", "TSS"], axis=1)

    Input_SNPs['start']=Input_SNPs['BP'].astype(int)-1
    Input_SNPs=Input_SNPs[['Chrom', 'start', 'BP']]

    #Annot_out=Input_SNPs['BP']   
    SNP_bed=BedTool.from_dataframe(Input_SNPs)
    if not os.path.exists(out_path+"chrX"):
        print("{} does not exist, creating".format(out_path+"chrX"))
        os.mkdir(out_path+"chrX/")
    for bed_for_annot in os.listdir(bed_for_annot_dir):
            Annot_out=Input_SNPs['BP']
            start_time_chr=time.time()
            Annot_bed=BedTool(bed_for_annot_dir+bed_for_annot)
            Annot_bed=Annot_bed.merge()
            Sample_name=bed_for_annot.split(".")[0]

            SNP_intersect=SNP_bed.intersect(Annot_bed)
            bp=[x.start + 1 for x in SNP_intersect]
            df_int=pd.DataFrame({'BP': bp, bed_for_annot: int(1)})
            Annot_out=pd.merge(Annot_out, df_int, how="left", on='BP')
            Annot_out.fillna(0, inplace=True)

            Annot_out=Annot_out.astype(int)
            Annot_out.to_csv(out_path+"chrX"+"/"+Sample_name, sep=" ", index=False)
            print("finished chrX"+" for "+Sample_name+" in "+str(time.time()-start_time_chr)+" seconds!")
    
    link_file=pd.DataFrame(columns=['Index', 'Annotation', 'Celltype', 'Tissue', 'Type', 'Category'])
    count=0

    for bed_for_annot in os.listdir(bed_for_annot_dir):
            #Annot_out=Input_SNPs['BP']
            start_time=time.time()

            Sample_name=bed_for_annot.split(".")[0]
            if(len(test_epimap(Sample_name, epimap_meta_path)) > 0):

                new_entry=pd.DataFrame({'Index': count, 'Annotation': Sample_name, 'Celltype':add_GROUP(Sample_name, epimap_meta_path), 
                'Tissue':add_infoline(Sample_name, epimap_meta_path), 'Type':"Enhancer", "Category":"Peaks"})

            else:
                new_entry=pd.DataFrame({'Index': count, 'Annotation': Sample_name, 'Celltype':"NA", 
                'Tissue':"NA", 'Type':"Enhancer", "Category":"Epimap"}, index=[count])
            link_file=pd.concat([link_file, new_entry])
            count=count+1
            print("finished link file for "+Sample_name+" in "+str(time.time()-start_time)+" seconds!")
    link_file['Celltype']=link_file.apply(lambda x: re.sub(" ", "_", x['Celltype']), axis=1) 
    link_file['Tissue']=link_file.apply(lambda x: re.sub(" ", "_", x['Tissue']), axis=1) 
    link_file.to_csv(out_path+"link_file.txt", sep=" ", index=False, header=['Index', 'Annotation', 'Celltype', 'Tissue', 'Type', 'Category'])



# def Garfield_annot_UK10K_line(input_SNP_path, out_path, bed_for_annot_dir, input_chr):
#     start_time=time.time()
#     with open(input_SNP_path+"chr"+str(input_chr), 'r') as f:
#         for line in f:
#             line=line.rstrip().split(" ")
#             bp=line[0];
#             current_chr="chr"+str(input_chr)
#             input_start=int(bp)-1
#             input_SNP=BedTool(str(current_chr)+" "+str(input_start)+" "+str(bp), from_string=True);

#             annot_out=[]

#             for bed_for_annot in os.listdir(bed_for_annot_dir):
#                 Annot_bed=BedTool(bed_for_annot_dir+bed_for_annot)
#                 SNP_intersect=input_SNP.intersect(Annot_bed)
#                 annot_out.append(str(len(SNP_intersect)))
            
#             with open(out_path+"chr"+str(input_chr), 'a') as out_file:
#                 out_file.write(str(bp)+" "+"".join(annot_out)+"\n")

#     print("finished chr"+str(input_chr)+" in "+str(time.time()-start_time)+" seconds!")

def Garfield_annot_UK10K_line(input_SNP_path, out_path, bed_for_annot, input_chr):
    start_time=time.time()
    with open(input_SNP_path+"chr"+str(input_chr), 'r') as f:
        for line in f:
            line=line.rstrip().split(" ")
            bp=line[0];
            current_chr="chr"+str(input_chr)
            input_start=int(bp)-1
            input_SNP=BedTool(str(current_chr)+" "+str(input_start)+" "+str(bp), from_string=True);

            

            Annot_bed=BedTool(bed_for_annot)
            SNP_intersect=input_SNP.intersect(Annot_bed)
            if len(SNP_intersect)==0 :
                annot_out=0
            else:
                annot_out=1
            with open(out_path, 'a') as out_file:
                out_file.write(str(bp)+" "+str(annot_out)+"\n")

    print("finished chr"+str(input_chr)+" in "+str(time.time()-start_time)+" seconds!")

#abc=Garfield_annot_UK10K_line("/Users/xg1/Downloads/garfield-data/maftssd/", 
#    "/Users/xg1/Downloads/garfield-data/Comparison_annotations/",
 #   "/Users/xg1/Downloads/ldsc/Comparison_enhancers/")

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('--Input_SNP_dir', metavar="<str>", type=str, required=True)
    parser.add_argument('--Input_bed', metavar="<str>", type=str, required=True)
    parser.add_argument('--Output_dir', metavar="<str>", type=str, required=True)
    
    args = parser.parse_args()
    return args

def main():
    args=parse_args()
    #Garfield_annot_UK10K_line(input_SNP_path=args.Input_SNP_dir, out_path=args.Output_dir, bed_for_annot=args.Input_bed, input_chr=args.Input_SNP_chr)
    Garfield_annot_UK10K(input_SNP_path=args.Input_SNP_dir, out_path=args.Output_dir, bed_for_annot_dir=args.Input_bed)
if __name__ == '__main__':

    main()
