#!/bin/bash
set -e

# get data
wget -Nc "https://i12g-gagneurweb.informatik.tu-muenchen.de/public/paper/drop_analysis/resource.tar.gz"
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
