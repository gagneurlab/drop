#!/bin/bash
set -e

# get data
resource_url="https://www.cmm.in.tum.de/public/paper/drop_analysis/resource.tar.gz"
wget -nc $resource_url
if [ ! -d Data ]; then
  rm -rf resource
	tar -zxvf resource.tar.gz
	mv resource Data
else
    echo "Data directory already exists, is not updated"
fi

# prepare data
cd Data
cp config_relative_wb1.8.yaml ../config.yaml
python fix_sample_anno.py

# unzip fasta
if [ ! -f "chr21.fa" ]; then gunzip chr21.fa.gz; fi
