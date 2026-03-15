#!/bin/bash


# Requires a server with 1Tb RAM!!!
# Use max 120 CPUs / threads in parallel
num_cores=120

# Define the branch that needs to be tested
branch_name=dev

# create test env first
env_path=./drop_test_env
mamba create -p ${env_path} "drop>=1.5.0" "pandoc>2.4" "r-base>=4.5"  
conda activate ${env_path}

# install new version from git branch 
pip install git+https://github.com/gagneurlab/drop@${branch_name}

# create fresh folder to run drop
run_folder=run_drop_${RANDOM}
mkdir -p ${run_folder}
cd ${run_folder}

# init drop and config
drop init
sed -re 's|^(root:).*|\1 "./"|' -i config.yaml
sed -re 's|^(htmlOutputPath:).*|\1 "./html"|' -i config.yaml
sed -re 's|^(sampleAnnotation:).*|\1 "/s/project/drop-analysis/kremer_checks/sample_annotation.tsv"|' -i config.yaml
sed -re 's|^(    v29:).*|\1 "/s/genomes/human/hg19/gencode29/gencode.v29lift37.sorted.gtf.gz"|' -i config.yaml
sed -re 's|^(genome:).*|\1 \n    ucsc: "/s/genomes/human/hg19/fasta/hg19.fa"|' -i config.yaml
sed -re 's|^(    qcVcf:).*|\1 "/s/public_webshare/public/paper/drop_analysis/resource/qc_vcf_1000G_hg19.vcf.gz"|' -i config.yaml

# run drop locally for all samples and modules
snakemake -j ${num_cores} -n \
    --rerun-triggers mtime \
    --default-resources ntasks=1 mem_mb=5000 gpu=0 \
    --jobs $num_cores

# run final checks on the results
Rscript ../check_kremer_full_drop.R


