# Detection of RNA Outlier Pipeline
[![Pipeline status](https://travis-ci.org/gagneurlab/drop.svg?branch=master)](https://travis-ci.org/gagneurlab/drop)
[![Version](https://img.shields.io/badge/Version-0.9.0-green.svg)](https://github.com/gagneurlab/drop/master)

The manuscript main file, supplementary figures and table can be found in the manuscript folder or in [protocol exchange](https://protocolexchange.researchsquare.com/article/993ff4a5-38ce-4261-902a-600dbd528ba2/v1).


[DROP documentation](https://drop-rna.readthedocs.io/en/latest/index.html)

<img src="drop_sticker.png" alt="drop logo" width="200" class="center"/>

## Dependencies
Programming languages:

+ python >= 3.6.7
     + pip >= 19.1
     + we recommend using a virtual environment e.g. anaconda
+ R >= 3.5 (https://www.r-project.org/). We recomment not to install R and its packages through conda.

### R packages
Bioconductor and base R packages need to be installed. The packages are listed in `drop/requirementsR.txt`. A script for installing these packages is provided. From the repository root just execute:
```
Rscript drop/installRPackages.R drop/requirementsR.txt
```
### Other packages
+ samtools >= 1.7 (https://www.htslib.org/download/)
+ bcftools (newest) (https://github.com/samtools/bcftools)
+ tabix (https://www.htslib.org/download/)
+ GATK (https://software.broadinstitute.org/gatk/)
+ graphviz (https://www.graphviz.org/)
+ pandoc (https://pandoc.org/)

## Installation
Make sure that all of the above listed [dependencies](#dependencies) are installed.
Then install DROP from github using `pip`. For this you need to recursively clone the repository with all its submodules first.
```
git clone https://github.com/gagneurlab/drop.git --recurse-submodules
```
Install DROP (activate your python environment if you are using one)
```
# conda activate drop_env
cd drop
pip install .
```
Alternatively, you can also install it directly without cloning
```
pip install git+https://github.com/gagneurlab/drop.git
```
Installation time (including all dependencies): ~ 1h

### Initialize a project
DROP projects are initialized in a separate directory dedicated to the analysis project. Calling the initialization command creates the necessary files.
```
cd <project/path>
drop init
``` 

## Set up the demo project

First, install the drop module according to [installation](#installation). 
Initialize the demo directory with a custom test project path.
In the following we will use `$HOME/drop_demo` as <project/path>.
```
cd $HOME/drop_demo
drop demo
```
This command downloads the demo data (10 BAM files and a VCF file with 10 samples from the GEUVADIS project subset to chromosome 21), initializes DROP and adapts the config file paths to your current project directory. 
Now the pipeline is ready to be executed using `snakemake`.
```
snakemake -n # dryrun
snakemake
```
Once the pipeline has run through, you will find the output in the `$HOME/drop_demo/Output`. It will consist of processed data, results and HTML pages. In order to view the complete HTML summary, open `$HOME/drop_demo/Output/htmlOutput/drop_demo_index.html` in the browser. Please note that this is a very small subset of the data used in the manuscript (100 samples and all chromosomes).

Expected runtime: 25 min

## Set up a custom project
Install the drop module according to [installation](#installation) and initialize the project in a custom project directory.
### Prepare the input data
Create a sample annotation that contains the sample IDs, file locations and other information necessary for the pipeline.
Edit the config file to set the correct file path of sample annotation and locations of non-sample specific input files. For these steps, please refer to the [documentation](https://drop-rna.readthedocs.io/en/latest/prepare.html).

### Execute the pipeline
Once these files are set up, you can execute a dry run from your project directory
```
snakemake -n
```
This shows you the rules of all subworkflows. Omit `-n` if you are sure that you want you execute all printed rules. You can also invoke single workflows explicitly e.g. for aberrant splicing with 
```
snakemake aberrantExpression -n
```

