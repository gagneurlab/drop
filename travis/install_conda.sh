#!/bin/bash
set -e

CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

if [ -e $HOME/miniconda/etc/profile.d/conda.sh ]
then
    echo "using cached miniconda"
else
    wget $CONDA_URL -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    source $HOME/miniconda/etc/profile.d/conda.sh
    hash -r
    conda config --set always_yes yes --set changeps1 no
fi

