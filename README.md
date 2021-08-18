# SNP_enrich_preprocess
Preprocessing for SNP tissue enrichment analysis methods

# Install dependencies into isolated environment

conda env create -n SNP_enrich --file environment.yaml

conda activate SNP_enrich

# Dependencies:
Requires bedtools2 for overlapping: 

$ wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz

$ tar -zxvf bedtools-2.29.1.tar.gz

$ cd bedtools2

$ make

On VM may require installing libraries:

sudo apt-get install build-essential

sudo apt-get install libz-dev

sudo apt install libbz2-dev

sudo apt install libclang-dev

sudo apt-get install liblzma-dev

sudo apt-get install bzip2
