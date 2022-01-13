#!/bin/bash

# helper file for RVC snakemake rule changeHeader

# 1 {input.bam}
# 2 {input.bai}
# 3 {params.ref}
# 4 {params.known_sites}
# 5 {params.ucsc2ncbi}
# 6 {params.ncbi2ucsc}
# 7 {log}
# 8 {resources.tmpdir}
# 9 {output.bqsr_table}

input_vcf=$1
repeat_mask=$2
ref=$3
tmpdir=$4
log=$5
output_vcf=$6

tmp_bed=$(mktemp)
tmp_bed="${tmp_bed}.bed"

vcf_chr=$(zgrep -v "^#" $input_vcf | cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" ||true )
bed_chr=$(cat $repeat_mask | cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" ||true  )

if [ $vcf_chr -eq 0  ] && [ $bed_chr -ne 0 ] #vcf has no chr, bed has chr
then
    # remove "chr" from the bed file
    sed -e "s/^chr//" $repeat_mask > $tmp_bed
elif [ $bed_chr -eq 0  ] && [ $vcf_chr -ne 0 ] #vcf has chr, bed has no chr
then
    # add "chr" to the bed file
    sed -e "s/^/chr/" $repeat_mask > $tmp_bed
else
    cp $repeat_mask $tmp_bed
fi

gatk IndexFeatureFile -I $tmp_bed
gatk --java-options -Djava.io.tmpdir=${tmpdir} VariantFiltration -R $ref \
-V $input_vcf --mask $tmp_bed -O $output_vcf 2> $log
