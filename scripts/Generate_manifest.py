## Generate manifest file for sumstats inputs:

import os
import sys
import json
import argparse
import gzip
import yaml
import pandas as pd
def main():
    # Args
    args = parse_args()
    Existing_enrichments_file=args.in_existing
    Input_sumstats_file=args.in_sumstats
    out_todo="format_commands.txt"
    out_done="finished_commands.txt"

    Finished_studies=[]
    with open(Existing_enrichments_file, 'r') as Existing_enrichments:
        for line in Existing_enrichments:
            fields=line.rstrip().split("/")
            #study_name=grep("parquet")
            study_handle=list(filter(lambda x:'results' in x, fields))
            if(study_handle ==[]):
                pass
            else:
                study_name=study_handle[0].split(".")[0]
                Finished_studies.append(study_name)

    todo_h=open(out_todo, "wb")
    done_h=open(out_done, "wb")

    with open(Input_sumstats_file, 'r') as Input_sumstats:
        for line in Input_sumstats:
            fields=line.rstrip().split("/")
            #study_name=grep("parquet")
            study_handle=list(filter(lambda x:'parquet' in x, fields))

            study_name=study_handle[0].split(".")[0]
            cmd=['python', "scripts/LDSC_format_single_sumstat.py", "--in_config", "configs/test_run_config.txt", "--in_study", study_name]
            cmd_str=" ".join([str(arg) for arg in cmd])
            print(cmd_str)
            if(study_name in Finished_studies):
                done_h.write((cmd_str + '\n').encode())
            else:
                todo_h.write((cmd_str + '\n').encode())

    ## python scripts/LDSC_format_single_sumstat.py --in_config configs/test_run_config.txt --in_study GCST002245


def parse_args():
    ''' Load command line args
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_existing', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_sumstats', metavar="<str>", type=str, required=True)


    args = parser.parse_args()
    return args


if __name__ == '__main__':

    main()