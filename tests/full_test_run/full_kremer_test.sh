#!/bin/bash
#
# Script to run a full fletched analysis on our wellknown Kremer et. al. data
# This ensures, that we can reidentify all known cases after each update.
#
# This full test should be run at every release and should not throw an error.
#
script_dir=$(realpath $(dirname $0))

# Requires a server with 1Tb RAM!!!
# Use max 200 CPUs / threads in parallel
num_cores=${1:-250}
branch_name=${2:-dev}
test_dir_root=${3:-$TMP}

# create fresh folder to run drop
run_folder=${test_dir_root}/drop-test/run_drop_${RANDOM}
mkdir -p ${run_folder}
cd ${run_folder}

# create test env first
env_path=./env/drop_test_env
mkdir -p env
mamba create -y -p ${env_path} \
    "drop>=1.5" "pandoc>=2.4" "r-base>=4.4" \
    "python>=3.12" "pip>=26.0" "yq>=3.4" \
    "pandas>=3.0"
conda activate ${env_path}

# install version from git branch
pip install --no-deps --force-reinstall \
    "git+https://github.com/gagneurlab/drop@${branch_name}"

# init drop and config
drop init

# adapt config to local setup
yq -yi '
  .root = "" |
  .htmlOutputPath = "html" |
  .sampleAnnotation = "/s/project/drop-analysis/kremer_checks/sample_annotation_80.tsv" |
  .geneAnnotation.v29 = "/s/genomes/human/hg19/gencode29/gencode.v29lift37.sorted.gtf.gz" |
  .genome.ucsc = "/s/genomes/human/hg19/fasta/hg19.fa" |
  .mae.groups = ["mae"] |
  .mae.qcVcf = "/s/public_webshare/public/paper/drop_analysis/resource/qc_vcf_1000G_hg19.vcf.gz" |
  .rnaVariantCalling = {} |
  .rnaVariantCalling.run = false
  ' config.yaml

# run drop locally for all samples and modules
snakemake -j 1 -n
snakemake --jobs ${num_cores} \
    --rerun-triggers mtime \
    --keep-going \
    --retries 3

# Run final check if we can recall all our findings
Rscript ${script_dir}/kremer_test.R

