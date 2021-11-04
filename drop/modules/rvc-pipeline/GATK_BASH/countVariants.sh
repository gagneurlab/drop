#!/bin/bash

# 1 {input.sample_vcf}
# 2 {wildcards.sample}
# 3 {output.vcf_counts}

sample_vcf=$1
sample_id=$2
output_counts=$3

total_variants=$(zgrep -vc "#" $sample_vcf)
pass_variants=$(zgrep -v "#" $sample_vcf |grep -c -P "PASS|\tMask\t")
pass_mask_variants=$(zgrep -v "#" $sample_vcf |grep -wc "PASS")

hom_total_variants=$(($(zgrep -v "#" $sample_vcf | grep -wc -F "1/1") +
$(zgrep -v "#" $sample_vcf | grep -wc -F "1|1")))

hom_pass_variants=$(($(zgrep -v "#" $sample_vcf | grep -w -F "1/1" |grep -c -P "PASS|\tMask\t" ) +
$(zgrep -v "#" $sample_vcf | grep -w -F "1|1" |grep -c -P "PASS|\tMask\t" )))

hom_pass_mask_variants=$(($(zgrep -v "#" $sample_vcf | grep -w -F "1/1" |grep -wc "PASS") +
$(zgrep -v "#" $sample_vcf | grep -w -F "1|1" |grep -wc "PASS")))

echo "Sample,GT,variants called,variants passing filter,variants passing filter and repeat mask" > $output_counts

echo "${sample_id},0/1,$(($total_variants - $hom_total_variants)),\
$(($pass_variants - $hom_pass_variants)),\
$(($pass_mask_variants - $hom_pass_mask_variants))" >> $output_counts

echo "${sample_id},1/1,$hom_total_variants,$hom_pass_variants,$hom_pass_mask_variants" >> $output_counts

