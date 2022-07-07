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

# CHEERS preprocessing

Generate epigenetic input:
	For EPIMAP: 

Generate consensus peak files:

```bash
$ python scripts/CHEERS_preprocessing/Create_full_epimap_consensus.py
```
	Which creates the consensus peaks in : 

```bash
$		../../tmp/Master_enhancers.sorted.merged.bed
```

Generate signals:

```bash
$	 python Generate_single_signal.py --Sample "BSS01668" --Peaks "../../tmp/Master_enhancers.sorted.merged.bed" --outdir "../../tmp/H3K27ac"
```

The output directory of the signals in consensus peak files can be passed to the CHEERS normalisation script.
		
	For BLUEPRINT: 
Generate consensus peak files:
```bash
$		python scripts/CHEERS_preprocessing/Generate_consensus_peaks.py --prefix BLUEPRINT --Peaks configs/BLUEPRINT_peaks.tsv --outdir ~/BLUEPRINT_peaks/
```

Generate signals in consensus peaks (for one sample):
(URLs are in configs/BLUEPRINT_signals.tsv)

```bash
$		python scripts/CHEERS_preprocessing/Generate_signal_from_BW_URL.py --Sample 0 --Peaks ~/BLUEPRINT_peaks/BLUEPRINT_Consensus_peaks.bed --BW_URL http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/bone_marrow/BM030613/band_form_neutrophil/ChIP-Seq/NCMLS/S00JGXH1.ERX651407.H3K27ac.bwa.GRCh38.20150529.bw --outdir /home/xg1/BLUEPRINT_peaks/ReadsInPeaks
```

Can be piped to parallel:
```bash
$	cat configs/BP_generate_signal_manifest.txt | parallel -j $NCORES
```

The output directory of the signals in consensus peak files can be passed to the CHEERS normalisation script.

#	Generating credible SNP set inputs:

We used the latest Open Targets finemapping results, the latest publicly available credible SNP sets can be downloaded at:
```bash
$ gs://open-targets-genetics-releases/22.02.01/v2d_credset
```

I’ve saved the credible SNPs as ~/Credible_SNP_sets/finemapping_220401.parquet/

Generate credible SNP sets in hg19 for EPIMAP:
```bash
$ python Liftover_hg38_to_hg19_single_study.py 
$   --Study_ID GCSTxxxx 
$   --input_credset ~/Credible_SNP_sets/finemapping_220401.parquet/
$   --Enrichment_outdir ~/CHEERS/Results/
$   --input_peak Normalised_signals.txt
$   --outdir ~/Credible_SNP_sets/Formatted_hg19/
```

This will also create a manifest file: ~/compute_CHEERS_enrichments.txt

Which can be used to run CHEERS compute enrichments in parallel:
```bash
$	cat ~/compute_CHEERS_enrichments.txt | parallel -j $NCORES
```

Generate credible SNP sets in hg38 for BLUEPRINT:
```bash
$ python No_liftover_single_study.py  --Study_ID GCST006979 --input_credset ~/Credible_SNP_sets/finemapping_220401.parquet/ --Enrichment_outdir /home/xg1/BLUEPRINT_peaks/Results/ --input_peak /home/xg1/BLUEPRINT_peaks/Normalised/BLUEPRINT_counts_normToMax_quantileNorm_euclideanNorm.txt --outdir ~/Credible_SNP_sets/Formatted_hg38/
```

This will also create a manifest file: ~/compute_CHEERS_enrichments.txt

Which can be used to run CHEERS compute enrichments in parallel:
```bash
$	cat ~/compute_CHEERS_enrichments.txt | parallel -j $NCORES
```
