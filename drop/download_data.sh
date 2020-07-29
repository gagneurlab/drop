#!/bin/bash
set -e

# get data
resource_url="https://www.cmm.in.tum.de/public/paper/drop_analysis/resource.tar.gz"
wget -nc $resource_url
if [ ! -d Data ]; then
	tar -zxvf resource.tar.gz
	mv resource Data
fi

# prepare data
cd Data
cp config_relative.yaml ../config.yaml
python fix_sample_anno.py

# unzip fasta
if [ ! -f "chr21.fa" ]; then gunzip chr21.fa.gz; fi
