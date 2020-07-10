#!/bin/bash
set -e

CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh"

if [ -e $HOME/miniconda/etc/profile.d/conda.sh ]
then
    echo "using cached miniconda"
else
    rm -rf $HOME/miniconda
    wget $CONDA_URL -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p $HOME/miniconda
    source $HOME/miniconda/etc/profile.d/conda.sh
    hash -r
    conda config --set always_yes yes --set changeps1 no
fi

