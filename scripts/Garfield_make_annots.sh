#!/bin/bash


INPUT_ANNOTATIONS=$1
OUTPUT_ANNOTATIONS=$2
DATA_DIR=$3


## Create annotation inputs for garfield:




# Reformat and concat into garfield format:

for i in {1..22} "X"; do 
    awk '{ print $1 }' ~/Downloads/garfield-data/maftssd/chr${i} > tmp.chr${i};
      for FILE in ~/Downloads/garfield-data/Epimap_EP_annots_new/chr${i}/* ; do 
        paste tmp.chr${i} <(awk '{ print $2}' <(tail -n +2 ${FILE})) > tmp0.chr${i} ;
        cp tmp0.chr${i} tmp.chr${i} ;
        done; 
    done

for f in {1..22} "X"; do 
  paste -d" " <(awk '{print $1}' tmp.chr${f} | sed 1d) <(awk '{$1=""; print $0}' tmp.chr${f} | awk '{ gsub("\t",""); print;}' | awk '{ gsub(" ",""); print;}'| sed 1d) > chr${f}; 
done

for f in {1..22} "X"; do 
  tail -n +2 chr${f}/BSS00196 > tmp0.chr${i} ;
  #cp /Users/xg1/Downloads/garfield-data/Epimap_EP_annots/chr${f}/BSS00196 /Users/xg1/Downloads/garfield-data/Epimap_EP_annots/chr${f}/tmp.chr${f};
  paste -d" " <(awk '{print $1}' tmp0.chr${i}   | sed 1d) <(awk '{$1=""; print $0}' tmp0.chr${i}  | awk '{ gsub("\t",""); print;}' | awk '{ gsub(" ",""); print;}'| sed 1d) > chr${f}.full; 
done


for FILE in ~/Downloads/garfield-data/test_concat/* ; do paste tmp.chrX <(awk '{ print $2}' <(tail -n +2 ${FILE})) > tmp0.chrX ; cp tmp0.chrX tmp.chrX ; done

paste -d" " <(awk '{print $1}' tmp.chrX | sed 1d) <(awk '{$1=""; print $0}' tmp.chrX | awk '{ gsub("\t",""); print;}' | awk '{ gsub(" ",""); print;}'| sed 1d) > chrX.test
