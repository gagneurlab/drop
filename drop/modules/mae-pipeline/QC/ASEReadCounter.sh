#!/bin/bash

# 1 {input.ncbi2ucsc}
# 2 {input.ucsc2ncbi}
# 3 {input.vcf_file_ucsc}
# 4 {input.vcf_file_ncbi}
# 5 {input.bam_file}
# 6 {wildcards.rna}
# 7 {input.fasta}
# 8 {config[mae][gatkIgnoreHeaderCheck]}
# 9 {output.counted}
# 10 {params.bcftools}
# 11 {params.samtools}
# 12 {params.gatk}
# 13 {input.script_mae}

ncbi2ucsc=$1
ucsc2ncbi=$2
vcf_file_ucsc=$3
vcf_file_ncbi=$4
bam_file=$5
rna_id=$6
fasta=$7
sanity=$8
output=$9
bcftools=${10}
samtools=${11}
gatk=${12}
script=${13}

# get number of UCSC chromosomes in BAM
bam_chr=$($samtools idxstats ${bam_file} | grep chr | wc -l)
if [ ${bam_chr} -ne 0 ]
then
    echo "use UCSC format"
    vcf_file=$vcf_file_ucsc
else
    echo "use NCBI format"
    vcf_file=$vcf_file_ncbi
fi

bash $script $ncbi2ucsc $ucsc2ncbi $vcf_file $bam_file \
    "qc--${rna_id}" $fasta $sanity $output $bcftools $samtools $gatk
