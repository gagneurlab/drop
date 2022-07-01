#!/bin/bash
set -e

# get data
resource_url="https://github.com/nickhsmith/drop_demo_data/archive/refs/heads/main.zip"
tmpdir="$(dirname "$(mktemp)")"
wget -nc -P $tmpdir $resource_url

# if the directory Data does not exist
if [ ! -d "Data" ]; then
	unzip "$tmpdir/main.zip" 
	mv drop_demo_data-main/Data Data
	rm -rf drop_demo_data-main
else
    echo "Data directory not empty, is not updated"
fi

# unzip fasta
cd ./Data
if [ ! -f "chr21.fa" ]; then gunzip chr21.fa.gz; fi
