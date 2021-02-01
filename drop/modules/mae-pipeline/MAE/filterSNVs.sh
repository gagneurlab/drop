#!/bin/bash
set -e

# 1 {input.ncbi2ucsc}
# 2 {input.ucsc2ncbi}
# 3 {input.vcf_file}
# 4 {wildcards.vcf}
# 5 {input.bam_file}
# 6 {output.snvs_filename}
# 7 {config[tools][bcftoolsCmd]}
# 8 {config[tools][samtoolsCmd]}

ncbi2ucsc=$1
ucsc2ncbi=$2
vcf_file=$3
vcf_id=$4
bam_file=$5
output=$6
bcftools=$7
samtools=$8

tmp=$(mktemp)

# $bcftools annotate --force -x INFO ${vcf_file} |\
$bcftools view  $vcf_file | grep -vP "^##INFO=" | awk -F'\t' 'BEGIN {OFS = FS} { if($1 ~ /^[^#]/){ $8 = "." }; print $0 }' |\
    $bcftools view -s ${vcf_id} -m2 -M2 -v snps -O z -o $tmp
$bcftools index -t $tmp

# compare and correct chromosome format mismatch
bam_chr=$($samtools idxstats ${bam_file} | cut -f1 | grep "^chr" | wc -l)
vcf_chr=$($bcftools index --stats $tmp   | cut -f1 | grep "^chr" | wc -l)

if [ ${vcf_chr} -eq 0  ] && [ ${bam_chr} -ne 0 ]
then
    echo "converting from NCBI to UCSC format"
    $bcftools annotate --rename-chrs $ncbi2ucsc $tmp | bgzip > ${output}
    rm ${tmp}
    rm ${tmp}.tbi
elif [ ${vcf_chr} -ne 0  ] && [ ${bam_chr} -eq 0 ]
then
    echo "converting from UCSC to NCBI format"
    $bcftools annotate --rename-chrs $ucsc2ncbi $tmp | bgzip > ${output}
    rm ${tmp}
    rm ${tmp}.tbi
else
    mv $tmp ${output}
    rm ${tmp}.tbi
fi

$bcftools index -t ${output}

