#!/bin/bash
set -e

# get data
resource_url="https://www.cmm.in.tum.de/public/paper/drop_analysis/resource_splice_merge.tar.gz"
tmpdir=$(dirname "$(mktemp -u)")
wget -nc -O "$tmpdir/resource.tar.gz" "$resource_url"

mkdir -p Data
if [ -z "$(ls Data)" ]; then
	tar -zxvf "$tmpdir/resource.tar.gz" -C .
	rm -rf Data
	mv resource Data
else
    echo "Data directory not empty, is not updated"
fi

# unzip fasta
cd ./Data
if [ ! -f "chr21.fa" ]; then gunzip chr21.fa.gz; fi
