#!/bin/bash


INPUT_SUMSTATS=$1
TRAIT_NAME=$2
DATA_DIR=$3

PVAL_DIR=${DATA_DIR}/pval/

## Reformat sumstats into garfield input:
## (Change lift_from and lift to accordingly)

python Garfield_make_pval.py --Input_sumstats ${INPUT_SUMSTATS} --Trait ${TRAIT_NAME}  --Output_dir ${PVAL_DIR}/${TRAIT_NAME} --Lift_from hg38 --Lift_to hg19

## Actually run Garfield:
