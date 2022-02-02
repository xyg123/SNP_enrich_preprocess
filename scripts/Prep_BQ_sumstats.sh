#!/usr/bin/env bash

Completed_enrichments=$1
## gs://genetics-portal-dev-analysis/xg1/LDSC_enrichments/*
Available_inputs=$2
## gs://genetics-portal-dev-analysis/xg1/Test_sumstat_inputs/*.parquet

NCORES=$3

gsutil -m ls -d $Completed_enrichments > existing_enrichments.txt
gsutil -m ls -d $Available_inputs > input_sumstats.txt

# writes the formatting commands to "format_commands.txt":

python scripts/Generate_manifest.py --in_existing existing_enrichments.txt --in_sumstats input_sumstats.txt



# pipe formatting sumstat commands to parallel:

cat format_commands.txt | parallel -j $NCORES --joblog logs/parallel.jobs.log


