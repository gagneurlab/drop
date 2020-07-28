# Detection of RNA Outlier Pipeline
[![Pipeline status](https://travis-ci.org/gagneurlab/drop.svg?branch=master)](https://travis-ci.org/gagneurlab/drop)
[![Version](https://img.shields.io/badge/Version-0.9.1-green.svg)](https://github.com/gagneurlab/drop/master)
[![Version](https://readthedocs.org/projects/gagneurlab-drop/badge/?version=latest)](https://gagneurlab-drop.readthedocs.io/en/latest)

The manuscript main file, supplementary figures and table can be found in the manuscript folder or in 
[protocol exchange](https://protocolexchange.researchsquare.com/article/993ff4a5-38ce-4261-902a-600dbd528ba2/v1).

<img src="drop_sticker.png" alt="drop logo" width="200" class="center"/>

## Installation
DROP is available on [bioconda](https://anaconda.org/bioconda/drop) for python 3.6 and above.
We recommend using a dedicated conda environment.

```
# create environment
conda create -n drop_env python=3.6
conda activate drop_env

# install drop
conda install -c bioconda drop
```
Installation time: ~ 10min

Test whether the pipeline runs through by setting up the demo dataset in an empty directory (e.g. ``~/drop_demo``).

```
mkdir ~/drop_demo
cd ~/drop_demo

# demo will download the necessary data and pipeline files
drop demo
```

The pipeline can be run using `snakemake` commands

```
snakemake -n # dryrun
snakemake
```

Expected runtime: 25 min

For more information on different installation options, check out the 
[documentation](https://gagneurlab-drop.readthedocs.io/en/latest/installation.html)

## Set up a custom project
Install the drop module according to [installation](#installation) and initialize the project in a custom project directory.
### Prepare the input data
Create a sample annotation that contains the sample IDs, file locations and other information necessary for the pipeline.
Edit the config file to set the correct file path of sample annotation and locations of non-sample specific input files.
The requirements are described in the [documentation](https://gagneurlab-drop.readthedocs.io/en/latest/prepare.html).

### Execute the pipeline
Once these files are set up, you can execute a dry run from your project directory
```
snakemake -n
```
This shows you the rules of all subworkflows. Omit `-n` if you are sure that you want you execute all printed rules. You can also invoke single workflows explicitly e.g. for aberrant splicing with 
```
snakemake aberrantExpression -n
```

## Datasets
The following publicly-available datasets of gene counts can be used as controls:

* 119 non-strand specific fibroblasts: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3887451.svg)](https://doi.org/10.5281/zenodo.3887451)

* 139 strand specific fibroblasts: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3963474.svg)](https://doi.org/10.5281/zenodo.3963474)

* 125 strand specific blood: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3963470.svg)](https://doi.org/10.5281/zenodo.3963470)

If you want to contribute with your own count matrices, please contact us: yepez at in.tum.de
