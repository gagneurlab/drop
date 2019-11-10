#!/bin/bash
set -e

cd $PROJECT_DIR
wget -nc "https://i12g-gagneurweb.informatik.tu-muenchen.de/public/paper/drop_analysis/resource.tar.gz"
tar -zxvf resource.tar.gz
mv resource Data
cd Data
pwd
ls
python fix_sample_anno.py
