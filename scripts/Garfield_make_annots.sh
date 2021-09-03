#!/bin/bash


INPUT_SUMSTATS=$1
TRAIT_NAME=$2
DATA_DIR=$3


## Create annotation inputs for garfield:




# Reformat and concat into garfield format:

for i in {1..22}; do awk '{ print $1 }' ~/Downloads/garfield-data/maftssd/chr${i} > tmp.chr${i};  for FILE in ~/Downloads/garfield-data/test/chr${i}/* ; do paste tmp.chr${i} <(awk '{ print $2}' <(tail -n +2 ${FILE})) > tmp0.chr${i} ; cp tmp0.chr${i} tmp.chr${i} ; done; done

for f in {1..22}; do cp  chr${f} chr${f}.tmp; paste -d" " <(awk '{print $1}' chr${f}.tmp | sed 1d) <(awk '{$1=""; print $0}' chr${f}.tmp | awk '{ gsub("\t",""); print;}' | awk '{ gsub(" ",""); print;}'| sed 1d) > chr${f}.test; done


 for FILE in ~/Downloads/garfield-data/test_concat/* ; do paste tmp.chrX <(awk '{ print $2}' <(tail -n +2 ${FILE})) > tmp0.chrX ; cp tmp0.chrX tmp.chrX ; done

paste -d" " <(awk '{print $1}' tmp.chrX | sed 1d) <(awk '{$1=""; print $0}' tmp.chrX | awk '{ gsub("\t",""); print;}' | awk '{ gsub(" ",""); print;}'| sed 1d) > chrX.test
