#!/bin/bash

# helper file for RVC snakemake rule indexReads

# 1 {input} sa.getFilePath(wildcards.sample)
# 2 {log}
# 3 {output.bam}
# 4 {output.bai}

input_bam=$1
log=$2
output_bam=$3
output_bai=$4

if [ ! -f "${input_bam}.bai" ]; then
    samtools index -b $input_bam 2> $log
fi

ln -s $input_bam $output_bam
ln -s "${input_bam}.bai" $output_bai
