#!/bin/bash
set -e

if [ ! -f $HOME/bcftools/bcftools ] || [ ! -d $HOME/htslib ]
then
    cd $HOME
    git clone git://github.com/samtools/htslib.git
    git clone git://github.com/samtools/bcftools.git
    cd bcftools
    make
fi

