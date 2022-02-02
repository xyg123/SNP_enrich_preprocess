#!/usr/bin/env bash


Completed_enrichments=$1
## gs://genetics-portal-dev-analysis/xg1/LDSC_enrichments
Available_inputs=$2
NCORES=$3

gsutil -m ls -d $Completed_enrichments > existing_enrichments.txt
gsutil -m ls -d $Available_inputs > input_sumstats.txt

python scripts/Generate_manifest.py --in_existing existing_enrichments.txt --in_sumstats input_sumstats.txt

cat format_commands.txt | parallel -j $NCORES --bar --joblog logs/parallel.jobs.log

#zcat Format_sumstats_commands | shuf | parallel -j $NCORES --joblog logs/parallel.jobs.log




# NCORES=$1

# set -euo pipefail

# python 3_make_commands.py | shuf | parallel -j $NCORES --bar --joblog logs/parallel.jobs.log

# echo COMPLETE

