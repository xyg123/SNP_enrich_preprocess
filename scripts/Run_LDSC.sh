#!/usr/bin/env bash
# Run LDSC 

Input_study=$1
LDSC_dir=$2
Annot_dir=$3

conda activate SNP_enrich
gsutil cp gs://genetics-portal-dev-analysis/xg1/rsid_sumstats/${Input_study} ../tmp/${Input_study}

conda activate ldsc

python munge_sumstats.py --sumstats ../tmp/${Input_study} --merge-alleles ${LDSC_dir}/w_hm3.snplist --out ${LDSC_dir}/munged_sumstats/${Input_study} --a1-inc
rm ../tmp/${Input_study}

#for FILE in ${Annot_dir}/*; do python ldsc.py --h2 ${LDSC_dir}/munged_sumstats/${Input_study} --ref-ld-chr ${LDSC_dir}/baselineLD/baselineLD.,${LDSC_dir}/Comparison_EP_concat/Epimap_EP. --frqfile-chr ${LDSC_dir}/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr ${LDSC_dir}/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --print-delete-vals --out ${LDSC_dir}/results/${Input_study}; done
for PART in {1..17}; do python ldsc.py --h2 ${LDSC_dir}/munged_sumstats/${Input_study} --ref-ld-chr ${LDSC_dir}/baselineLD/baselineLD.,${Annot_dir}/part_${PART}/Epimap. --frqfile-chr ${LDSC_dir}/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr ${LDSC_dir}/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --print-delete-vals --out ${LDSC_dir}/results/${Input_study}; done

# Run script to merge results:

