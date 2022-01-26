from google.cloud import bigquery
client = bigquery.Client()

#Input_study="GCST000568"
Input_study=["GCST005536", "GCST005413", "GCST004131", "GCST005569", "GCST002245"]
## LOAD Job: load from GCS to BQ (table_id)
    # Would it be possible to make this table temporary? Or delete itself automatically after 1 week?

Input_sumstats_path="gs://genetics-portal-dev-sumstats/unfiltered/gwas"

for study in Input_study:
    Input_study_URI=Input_sumstats_path+"/"+study+".parquet/*.parquet"
    table_id = "open-targets-genetics-dev.Jack_sumstats."+study

    load_job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.PARQUET,)
    load_job = client.load_table_from_uri(
        Input_study_URI, table_id, job_config=load_job_config
    )  # Make an API request.
    load_job.result()  # Waits for the job to complete.

    destination_table = client.get_table(table_id)  # Make an API request.
    print("Loaded {} rows.".format(destination_table.num_rows))

# Query Job

for study in Input_study:

    table_id = "open-targets-genetics-dev.Jack_sumstats."+study
    rsID_table="open-targets-genetics-dev.rsID_sumstats."+study

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

for study in Input_study:
    rsID_table="open-targets-genetics-dev.rsID_sumstats."+study

    rsID_GCS_bucket="gs://genetics-portal-dev-analysis/xg1/rsid_sumstats"
    rsID_GCS_URI=rsID_GCS_bucket+"/{0}.txt.gz".format(study)

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

