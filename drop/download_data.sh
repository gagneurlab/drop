#!/bin/bash
set -e

# get data
resource_url="https://github.com/nickhsmith/drop_demo_data/archive/refs/heads/main.zip"
tmpdir="$(dirname "$(mktemp)")"
wget -nc -P $tmpdir $resource_url
mkdir -p Data
if [ -z "$(ls Data)" ]; then
    unzip "$tmpdir/main.zip"
    rm -rf Data
    mv drop_demo_data-main/Data Data
    rm -rf drop_demo_data-main
else
    echo "Data directory not empty, is not updated"
fi

