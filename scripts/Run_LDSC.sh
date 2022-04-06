#!/usr/bin/env bash
# Run LDSC :
# IMPORTANT: conda activate ldsc 

Input_study=$1
LDSC_dir=$2
Annot_dir=$3

#conda activate SNP_enrich
mkdir ../tmp/${Input_study}
gsutil rsync gs://genetics-portal-dev-analysis/xg1/rsid_sumstats/${Input_study}/ ../tmp/${Input_study}/

#conda activate ldsc
cd ../tmp/${Input_study}

files=( * )
if [[ ${#files[@]} -gt 1 ]]; then
    mv ${files[0]} ${Input_study}.txt.gz
    max=$(expr ${#files[@]} - 1)
    for i in $(seq 1 1 $max); do 
	    zcat ${Input_study}.txt.gz <(zcat ${files[${i}]} | tail -n +2 | gzip ) | gzip > tmp.txt.gz ;
            mv tmp.txt.gz ${Input_study}.txt.gz ;
    done
else
    mv ${files[0]} ${Input_study}.txt.gz
fi

cd ~/SNP_enrich_preprocess

#mkdir ${LDSC_dir}/munged_sumstats/

python ${LDSC_dir}/munge_sumstats.py --sumstats ../tmp/${Input_study}/${Input_study}.txt.gz --merge-alleles ${LDSC_dir}/w_hm3.snplist --out ${LDSC_dir}/munged_sumstats/${Input_study} --a1-inc --chunksize 50000

rm ../tmp/${Input_study}/*

#for FILE in ${Annot_dir}/*; do python ldsc.py --h2 ${LDSC_dir}/munged_sumstats/${Input_study} --ref-ld-chr ${LDSC_dir}/baselineLD/baselineLD.,${LDSC_dir}/Comparison_EP_concat/Epimap_EP. --frqfile-chr ${LDSC_dir}/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr ${LDSC_dir}/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --print-delete-vals --out ${LDSC_dir}/results/${Input_study}; done

mkdir ${LDSC_dir}/results/${Input_study}
for PART in {1..4}; do python ${LDSC_dir}/ldsc.py --h2 ${LDSC_dir}/munged_sumstats/${Input_study}.sumstats.gz --ref-ld-chr ${LDSC_dir}/baselineLD/baselineLD.,${Annot_dir}/part_${PART}/Epimap.part_${PART}. --frqfile-chr ${LDSC_dir}/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr ${LDSC_dir}/weights_hm3_no_hla/weights. --overlap-annot --print-coefficients --print-delete-vals --out ${LDSC_dir}/results/${Input_study}/${Input_study}.part_${PART}; done

# Run script to merge results:

#bash Run_LDSC.sh GCST002245 ~/ldsc  ~/LDSC_full_epimap_parts
