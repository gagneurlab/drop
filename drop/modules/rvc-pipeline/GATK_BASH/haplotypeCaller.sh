#!/bin/bash

# helper file for RVC snakemake rule changeHeader

# 1 {input.bam}
# 2 {input.bai}
# 3 {params.ref}
# 4 {params.dbSNP}
# 5 {params.ucsc2ncbi}
# 6 {params.ncbi2ucsc}
# 7 {log}
# 8 {resources.tmpdir}
# 9 {output}
# 10 {hcArgs}

input_bam=$1
input_bai=$2
ref=$3
dbSNP=$4
ucsc2ncbi=$5
ncbi2ucsc=$6
log=$7
tmpdir=$8
output_gVCF=$9
hcArgs=${10}

# use samtools and bcftools to identify whether the bam file and
# the dbSNP file are in the same chr format.
# Use || true to avoid erros on an empty grep search

bam_chr=$(samtools idxstats $input_bam | cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" || true)
vcf_chr=$(bcftools index --stats $dbSNP| cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" || true)

if [ $bam_chr -eq 0  ] && [ $vcf_chr -ne 0 ] #bam has no chr, vcf has chr styling
then
    echo "converting dbSNP from NCBI to UCSC format" | tee $log
    tmp_vcf="$(mktemp).vcf.gz"
    bcftools annotate --rename-chrs $ucsc2ncbi $dbSNP | bgzip > "${tmp_vcf}"
    tabix -f "${tmp_vcf}"
elif [ $bam_chr -ne 0  ] && [ $vcf_chr -eq 0 ] #bam has chr, vcf has no chr styling
then
    echo "converting dbSNP from UCSC to NCBI format" | tee $log
    tmp_vcf="$(mktemp).vcf.gz"
    bcftools annotate --rename-chrs $ncbi2ucsc $dbSNP | bgzip > "${tmp_vcf}"
    tabix -f "${tmp_vcf}"
else
    echo "chromosome styles match" |tee $log
    tmp_vcf=$dbSNP

fi

echo "starting HaplotypeCaller"
# using the tmp known_sites vcf use HaplotypeCaller
gatk --java-options -Djava.io.tmpdir=${tmpdir} HaplotypeCaller -I $input_bam -R $ref \
--dont-use-soft-clipped-bases -stand-call-conf 20.0 --dbsnp "${tmp_vcf}" \
--output-mode EMIT_ALL_CONFIDENT_SITES -ERC GVCF $hcArgs  \
-O $output_gVCF 2>&1 | tee -a $log
