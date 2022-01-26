import os

##list_of_sumstats=os.system("bq ls --project_id=open-targets-ukbb --max_results=100000 outcomes")

# The outputs are wrote to stdout

# Test with a pre-defined list of studies:
sumstat_studies=["GCST000568", "GCST000569", "GCST000571", "GCST000612", "GCST000679"]
os.system()

# Loading data into BQ shouldn't be affected by this

# Neither should running queries?