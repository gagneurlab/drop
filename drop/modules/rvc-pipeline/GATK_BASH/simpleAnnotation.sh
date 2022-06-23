#!/bin/bash

# helper file for RVC snakemake rule changeHeader

# 1 {input.vcf}
# 2 {input.gtf}
# 3 {output.vcf}

input_vcf=$1
input_gtf=$2
output_vcf=$3

tmp_bed=$(mktemp)
tmp_bed="${tmp_bed}.bed"


if [[ $input_gtf == *.gz ]];
    then
        zcat $input_gtf | cut -f1,4,5,9 > $tmp_bed
    else
        cut -f1,4,5,9 $input_gtf > $tmp_bed
fi

# remove all ';' characters from the annotation to be parsed as 1 field
sed -i "s/;//g" $tmp_bed

#grep the first column (no headers) for chr1-9,chr10-19,chr20-22
vcf_chr=$(zgrep -v "^#" $input_vcf | cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$"||true )
bed_chr=$(cat $tmp_bed | cut -f1 |  grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$"||true )

if [ $vcf_chr -eq 0  ] && [ $bed_chr -ne 0 ] #vcf has no chr, bed has chr
then
    # remove "chr" from the bed file
    sed -i -e "s/^chr//" $tmp_bed
elif [ $bed_chr -eq 0  ] && [ $vcf_chr -ne 0 ] #vcf has chr, bed has no chr
then
    # add "chr" to the bed file
    sed -i -e "s/^/^chr/" $tmp_bed
fi

bcftools annotate -a $tmp_bed -c CHROM,FROM,TO,GENE \
 -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') $input_vcf | bgzip > $output_vcf

tabix $output_vcf
