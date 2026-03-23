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
env_path=${TMP}/env/env-$(basename ${run_folder})
mkdir -p $(dirname ${env_path})
wget "https://github.com/gagneurlab/drop/blob/${branch_name}/environment.yml"
yq -r '.dependencies[] | strings' environment.yml > environment.txt
mamba create -y -p ${env_path} --file ./environment.txt
source $(dirname ${CONDA_EXE})/activate ${env_path}

# Report locations in the end
report_folder_on_exit(){
    echo "The output data is in: ${run_folder} and the env is in ${env_path}"
}
trap report_folder_on_exit EXIT
trap report_folder_on_exit INT


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

