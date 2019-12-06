#!/bin/bash
set -e

PROJECT_DIR=$1
echo $PROJECT_DIR
mkdir -p $PROJECT_DIR
cd $PROJECT_DIR
wget -nc "https://i12g-gagneurweb.informatik.tu-muenchen.de/public/paper/drop_analysis/resource.tar.gz"
tar -zxvf resource.tar.gz
rm -rf Data
mv resource Data
cd Data
python fix_sample_anno.py
gunzip chr21.fa.gz
samtools faidx chr21.fa
gatk CreateSequenceDictionary --REFERENCE chr21.fa

