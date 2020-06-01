#!/bin/bash
set -e

# get data
resource_url="https://i12g-gagneurweb.in.tum.de/public/paper/drop_analysis/resource.tar.gz"
wget -Nc --no-check-certificate $resource_url
tar -zxvf resource.tar.gz
rm -rf Data
mv resource Data

# prepare data
cd Data
python fix_sample_anno.py
gunzip chr21.fa.gz
samtools faidx chr21.fa

# copy config
cp config_relative.yaml ../config.yaml
