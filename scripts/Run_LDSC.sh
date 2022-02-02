#!/usr/bin/env bash
# Run LDSC 

Input_study=$1
LDSC_dir=$2
Annot_dir=$3

#conda activate SNP_enrich
#gsutil cp gs://genetics-portal-dev-analysis/xg1/rsid_sumstats/${Input_study}.txt.gz tmp/${Input_study}.txt.gz

#source /home/xg1/miniconda3/bin/activate ldsc

python ${LDSC_dir}/munge_sumstats.py --sumstats tmp/${Input_study}.txt.gz --merge-alleles ${LDSC_dir}/w_hm3.snplist --out ${LDSC_dir}/munged_sumstats/${Input_study} --a1-inc
#rm ../tmp/${Input_study}.txt.gz

#for FILE in ${Annot_dir}/*; do python ldsc.py --h2 ${LDSC_dir}/munged_sumstats/${Input_study} --ref-ld-chr ${LDSC_dir}/baselineLD/baselineLD.,${LDSC_dir}/Comparison_EP_concat/Epimap_EP. --frqfile-chr ${LDSC_dir}/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr ${LDSC_dir}/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --print-delete-vals --out ${LDSC_dir}/results/${Input_study}; done
mkdir ${LDSC_dir}/results/${Input_study}
for PART in {1..17}; do python ${LDSC_dir}/ldsc.py --h2 ${LDSC_dir}/munged_sumstats/${Input_study}.sumstats.gz --ref-ld-chr ${LDSC_dir}/baselineLD/baselineLD.,${Annot_dir}/part_${PART}/Epimap.part_${PART}. --frqfile-chr ${LDSC_dir}/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr ${LDSC_dir}/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --print-delete-vals --out ${LDSC_dir}/results/${Input_study}/${Input_study}.part_${PART}; done

# Run script to merge results:

