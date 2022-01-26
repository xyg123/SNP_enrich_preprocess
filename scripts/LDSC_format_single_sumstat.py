from google.cloud import bigquery
import os
import sys
import json
import argparse
import gzip
import configparser
import pandas as pd


def main():
    # Load args
    args = parse_args()
    In_config=args.in_config
    Input_study=args.in_study

    Configs = configparser.ConfigParser()
    Configs.read(In_config)

    client = bigquery.Client()

    ## LOAD Job: load from GCS to BQ (table_id)
        # Would it be possible to make this table temporary? Or delete itself automatically after 1 week?

    Input_sumstats_path=Configs.get("config", "Input_sumstats_GCS")


    Input_study_URI=Input_sumstats_path+"/"+Input_study+".parquet/*.parquet"
    temp_BQ_sumstats=Configs.get("config", "Temp_BQ_sumstats")
    table_id = temp_BQ_sumstats+"."+Input_study
    print(table_id)
    load_job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.PARQUET,)
    load_job = client.load_table_from_uri(
        Input_study_URI, table_id, job_config=load_job_config
    )  # Make an API request.
    load_job.result()  # Waits for the job to complete.

    destination_table = client.get_table(table_id)  # Make an API request.
    print("Loaded {} rows.".format(destination_table.num_rows))

    # Query Job

    table_id = Configs.get("config", "Temp_BQ_sumstats")+"."+Input_study
    rsID_table = Configs.get("config", "RSID_BQ_sumstats")+"."+Input_study

    query_job_config = bigquery.QueryJobConfig(destination=rsID_table)
    query = """
        WITH SNP_info AS (
        SELECT
            CONCAT(CAST(chrom AS string), CAST(pos AS string), CAST(ref AS string), CAST(alt AS string)) AS identifier,
            ref,
            alt,
            n_total,
            pval,
            eaf,
            beta
        FROM
            `{0}` )
        SELECT
        rs_id AS RSID, ref AS A1, alt AS A2, n_total AS N, pval AS P, eaf AS EAF, beta AS BETA
        FROM
        SNP_info
        JOIN (
        SELECT
            CONCAT(CAST(chr_id AS string), CAST(position AS string), CAST(ref_allele AS string), CAST(alt_allele AS string)) AS identifier,
            rs_id
        FROM
            `open-targets-genetics.210608.variants` ) variants
        USING(identifier)
    """.format(table_id)

    query_job = client.query(query, job_config=query_job_config)
    query_job.result()

    # Extract Job

    rsID_GCS_bucket=Configs.get("config", "Formatted_sumstats_GCS")
    rsID_GCS_URI=rsID_GCS_bucket+"/{0}.txt.gz".format(Input_study)

    extract_job_config = bigquery.ExtractJobConfig() 
    extract_job_config.field_delimiter = '\t'
    extract_job_config.compression='GZIP'

    extract_job = client.extract_table(
        rsID_table,
        rsID_GCS_URI,
        # Location must match that of the source table.
        location="EU",
        job_config=extract_job_config
    )  # API request
    extract_job.result()  # Waits for job to complete.

    print(
        "Exported {} to {}".format(rsID_table, rsID_GCS_URI)
    )


def parse_args():
    ''' Load command line args
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_config', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_study', metavar="<str>", type=str, required=True, help=("Study ID of input sumstats"))

    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()