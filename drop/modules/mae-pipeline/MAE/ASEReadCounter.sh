#!/bin/bash

# 1 {input.ncbi2ucsc}
# 2 {input.ucsc2ncbi}
# 3 {input.vcf_file}
# 4 {input.bam_file}
# 5 {wildcards.vcf}--{wildcards.rna}
# 6 {input.fasta}
# 7 {config[mae][gatkIgnoreHeaderCheck]}
# 8 {output.counted}
# 9 {params.bcftools}
#10 {params.samtools}
#11 {params.gatk}

ncbi2ucsc=$1
ucsc2ncbi=$2
vcf_file=$3
bam_file=$4
mae_id=$5
fasta=$6
sanity=$7
output=$8
bcftools=$9
samtools=${10}
gatk=${11}

tmp=$(mktemp)
header="contig\tposition\tvariantID\trefAllele\taltAllele\t"
header+="refCount\taltCount\ttotalCount\tlowMAPQDepth\t"
header+="lowBaseQDepth\trawDepth\totherBases\timproperPairs"
echo -e $header >> $tmp

# get chr format
bam_chr=$($samtools idxstats ${bam_file} | grep -vP "\t0\t0" | cut -f1) # only chr with coverage
vcf_chr=$($bcftools view ${vcf_file} | cut -f1 | grep -v '#' | uniq)
if [ "$(echo ${vcf_chr} | grep -c 'chr')" -eq 0 ]; then
  echo "use NCBI format"
  canonical=$ncbi2ucsc
else
  echo "use UCSC format"
  canonical=$ucsc2ncbi
fi
# subset from canonical chromosomes
chr_subset=$(comm -12 <(cut -f1 -d" " ${canonical} | sort -u) <(echo ${vcf_chr} | xargs -n1 | sort -u))
chr_subset=$(comm -12 <(echo ${bam_chr} | xargs -n1) <(echo ${chr_subset} | xargs -n1) | uniq)

for chr in $chr_subset; do
  echo $chr
  $gatk ASEReadCounter \
    -R ${fasta} \
    -I ${bam_file} \
    -V ${vcf_file} \
    -L ${chr} \
    --verbosity ERROR \
    --QUIET true \
    --disable-sequence-dictionary-validation ${sanity} |
    tail -n+2 >>$tmp
done

echo $mae_id
cat $tmp | awk -v id="${mae_id}" \
  -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "ID"} NR>1{print $0, id}' |
  bgzip >${output}
rm ${tmp}

zcat ${output} | head
