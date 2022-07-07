import os
import sys
import argparse
import pandas as pd


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--All_signals', metavar="<str>", help=("All input signals"), type=str, required=True)
    #parser.add_argument('--Celltype', metavar="<str>", help=("Cell type to average across (NOT TISSUE)"), type=str, required=True)
    parser.add_argument('--Metadata', metavar="<str>", help=("Metadata to match input signal samples to cell type"), type=str, required=True)
    parser.add_argument('--Outdir', metavar="<str>", help=("Output directory"), type=str, required=True)
    #parser.add_argument('--input_peak', metavar="<str>", help=("Input signal at peaks file"), type=str, required=True)
    #parser.add_argument('--outdir', metavar="<str>", help=("Output for hg19 SNPs"), type=str, required=True)
    
    args = parser.parse_args()
    return args

def main():
    args=parse_args()

    
    # Read in All_signals
    print("Reading in full signals file...")
    Full_signals=pd.read_csv(args.All_signals, sep="\t")
    
    Samples=Full_signals.columns[3:]
    New_sample_names=list()
    for i in range(0, len(Samples)):
        New_sample_names.append(Samples[i].split("_")[0])
    New_sample_names=["chr", "start", "end"]+New_sample_names
    Full_signals.columns=New_sample_names
    #Full_signals
    
    
    # Identify which samples belong to Celltype

    Epimap_meta=pd.read_csv(args.Metadata, sep="\t")
    # Drop "Cancer" cell types
    Epimap_meta=Epimap_meta.loc[~Epimap_meta["GROUP"].str.contains("Cancer")]
    # Also drop "Other", it confuses the combination step.
    Epimap_meta=Epimap_meta.loc[~Epimap_meta["GROUP"].str.contains("Other")]

    # Also drop Samples that have been stimulated/perturbed:
    Epimap_meta=Epimap_meta.loc[Epimap_meta["perturb"].isnull()]

    Epimap_meta=Epimap_meta.loc[Epimap_meta["id"].isin(list(Full_signals.columns))]


    # Loop through each tissue and cell type:
    Tissue_types=Epimap_meta["GROUP"].unique()
    for i in range(0, len(Tissue_types)):
        Sub_celltypes=Epimap_meta.loc[Epimap_meta["GROUP"].str.contains(Tissue_types[i])]["infoline"].unique()
        Tissue_signals=Full_signals[["chr", "start", "end"]]
        print("Processing "+str(Tissue_types[i])+"...")
        for j in range(0, len(Sub_celltypes)):
            Matched_samples=Epimap_meta.loc[Epimap_meta["infoline"] == Sub_celltypes[j]]["id"]
            Matched_signals=Full_signals[list(Matched_samples)]
            if(len(Matched_samples) == 1):
                # Just add to Tissue out matrix
                Tissue_signals=pd.concat([Tissue_signals.reset_index(drop=True), Matched_signals.reset_index(drop=True)], axis=1)
            else:
                # Apply an average mean function across all rows
                Average_signals=Matched_signals.apply(lambda x: x.mean(), axis=1)
                # Add to Tissue out matrix
                Tissue_signals=pd.concat([Tissue_signals.reset_index(drop=True), Average_signals.reset_index(drop=True)], axis=1)
            Sub_celltypes[j]=Sub_celltypes[j].replace(" ", "_")
        Tissue_signals.columns=["chr", "start", "end"]+list(Sub_celltypes)
        # write to file: Outdir+Tissue_type+.txt
        Tissue_types[i]=Tissue_types[i].replace(" ", "_")
        Tissue_signals.to_csv(args.Outdir+str(Tissue_types[i])+".txt", sep="\t")

                
    


if __name__ == '__main__':
    
    main()

# python scripts/CHEERS_preprocessing/Average_signal_by_celltype.py \
#   --All_signals ~/Downloads/Full_Epimap_H3K27ac_testfile.txt \
#   --Metadata configs/main_metadata_table.tsv \
#   --Outdir ~/Downloads/Epimap_by_tissue/
