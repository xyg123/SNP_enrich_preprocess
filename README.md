# SNP_enrich_preprocess
Preprocessing for SNP tissue enrichment analysis methods

## Install dependencies into isolated environment
```bash
conda env create -n SNP_enrich --file environment.yaml

conda activate SNP_enrich
```
## Dependencies:
Requires bedtools2 in the $PATH for overlapping: 
```bash
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz

$ tar -zxvf bedtools-2.29.1.tar.gz

$ cd bedtools2

$ make
```
On VM may require installing libraries:

```bash
sudo apt-get install build-essential

sudo apt-get install libz-dev

sudo apt install libbz2-dev

sudo apt install libclang-dev

sudo apt-get install liblzma-dev

sudo apt-get install bzip2
```

# LDSC Preprocessing:
## Download all available Epimap enhancer/promoter files and merge them into one bed file:

```bash
$ python LDSC_get_all_epimap.py
```

In future, please check if bed files are hg19/Hg38 and liftover if needed. (Sumstats from GCS are Hg38)

## LDSC partitioned heritability:

Edit config file to specify the GCS buckets for input sumstats, temporary processing folders, and already existing results.

```bash
$ vim configs/configs.txt
```
### Formatting sumstats:
The following will process all summary stats which do not already have LDSC results:

```bash
$ NCORES = number of cores

$ scripts/Prep_BQ_sumstats.sh gs://genetics-portal-dev-analysis/xg1/LDSC_enrichments/* gs://genetics-portal-dev-analysis/xg1/Test_sumstat_inputs/*.parquet $NCORES
```

### Make commands to run LDSC:

```bash
$ INPUT_STUDY_PATH=gs://genetics-portal-dev-analysis/xg1/rsid_sumstats/

$ Make_LDSC_manifest.sh $INPUT_STUDY_PATH
```
### Run LDSC:
Commands are written to LDSC_studies_to_run.txt

Pipe to parallel:

```bash
$ cat LDSC_studies_to_run.txt | parallel -j $NCORES --joblog logs/parallel.jobs.log
```


