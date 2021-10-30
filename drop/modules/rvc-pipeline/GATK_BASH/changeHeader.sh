#!/bin/bash

# helper file for RVC snakemake rule changeHeader

# 1 {input.bam}
# 2 {input.bai}
# 3 {params.sample}
# 4 {log}
# 5 {output.bam}
# 6 {output.bai}
# 7 {output.newHeader}
input_bam=$1
input_bai=$2
sample=$3
log=$4
output_bam=$5
output_bai=$6
output_newHeader=$7


samtools view -H $input_bam |grep "^@RG" |grep SM |head -1|
while read header ; do
    for i in $header; do
        if [[ $i == "SM:"* ]]; then
            internalHeader=${i:3}
            break
        else
            internalHeader=""
        fi
    done

    if [[ $internalHeader == $sample ]]; then
        echo "Internal Header $internalHeader matches $sample" |tee $log
        echo "Internal Header is designated: $internalHeader"  |tee -a $log
        echo "SampleID is $sample" |tee -a $log
        ln -s  $input_bam $output_bam
        ln -s  $input_bai $output_bai
        echo "Done Linking files"
    else
        echo "WARNING"
        echo "Internal Header is designated: $internalHeader" |tee $log
        echo "SampleID is $sample"  |tee -a $log
        echo "Forcing $internalHeader to match $sample" |tee -a $log

        samtools view -H $input_bam > $output_newHeader
        echo $output_newHeader

        # sed using regEx in place substitiute 'SM:' followed by any thing that isn't tab or newLine. and then replace it with the sampleID and the delimiter (tab or newLine) that matched in the 1st expression.
        sed -E -i "s/(SM:[^\t|\n]*)(\t|\n*)/SM:${sample}\2/" ${output_newHeader}

        samtools reheader $output_newHeader $input_bam > $output_bam

        samtools index -b $output_bam

    fi
    echo "new header can be found here:$output_newHeader" |tee -a $log
done
