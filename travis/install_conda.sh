#!/bin/bash
set -e

cd $HOME
CONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"

if [ -d $MINICONDA_DIR ] && [ -e $CONDA_SCRIPT ]
then
    echo "using cached miniconda"
else
    rm -rf $MINICONDA_DIR
    wget $CONDA_URL -O miniconda.sh
    bash miniconda.sh -b -p $MINICONDA_DIR
    source $CONDA_SCRIPT
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda create -q -n drop python=$TRAVIS_PYTHON_VERSION
    conda env update -f conda.recipe/environment.yaml
    conda activate drop
    java -version
    gatk --help
    samtools --version
    bcftools --version
fi

