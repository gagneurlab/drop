# Detection of RNA Outlier Pipeline
[![DROP pipeline status](https://github.com/gagneurlab/drop/workflows/Build/badge.svg?branch=master)](https://github.com/gagneurlab/drop/actions?query=workflow%3ABuild)
[![Version](https://img.shields.io/github/v/release/gagneurlab/drop?include_prereleases)](https://github.com/gagneurlab/drop/releases)
[![Version](https://readthedocs.org/projects/gagneurlab-drop/badge/?version=latest)](https://gagneurlab-drop.readthedocs.io/en/latest)


The detection of RNA Outliers Pipeline (DROP) is an integrative workflow to detect aberrant expression, aberrant splicing, and mono-allelic expression from raw RNA sequencing data. 

The manuscript is available in [Nature Protocols](https://www.nature.com/articles/s41596-020-00462-5).

This [website](https://cmm.in.tum.de/public/paper/drop_analysis/webDir/drop_analysis_index.html) contains the different reports of the Geuvadis demo dataset described in the paper.

This [video](https://www.youtube.com/watch?v=XvgjiFQClhM&t=2761s) presents the tools used in DROP and their application to rare disease diagnostics.


<img src="drop_sticker.png" alt="drop logo" width="200" class="center"/>


## Quickstart
DROP is available on [bioconda](https://anaconda.org/bioconda/drop).
We recommend using a dedicated conda environment (`drop_env` in this example). Installation time: ~ 15 min.
```
mamba create -n drop_env -c conda-forge -c bioconda drop --override-channels
```


In the case troubles with mamba or conda, we recommend using the fixed `DROP_<version>.yaml` installation file we make available on our [public server](https://www.cmm.in.tum.de/public/paper/drop_analysis/). Install the current version and use the full path in the following command to install the conda environment `drop_env`
```
mamba env create -f DROP_1.4.0.yaml

```

Test installation with demo project
```
conda activate drop_env
mkdir ~/drop_demo
cd ~/drop_demo
drop demo
```

The pipeline can be run using [snakemake](https://snakemake.readthedocs.io/) commands
```
snakemake -n # dryrun
snakemake --cores 1
```

Expected runtime: 25 min

For more information on different installation options, refer to the
[documentation](https://gagneurlab-drop.readthedocs.io/en/latest/installation.html)

## What's new


Version 1.4.0 fixes some bugs regarding the split reads counting for FRASER, which only affects cohorts containing samples with different strands (stranded forward, stranded reverse or unstranded). If your cohort contained samples with different strands, please rerun the AS module using version 1.4.0. In addition, due to snakemake updates affecting wBuild and the way we installed FRASER, installing DROP 1.3.3 no longer works. Version 1.3.4 fixes the FRASER version to ensure reproducibility and fixes certain scripts affected by the snakemake update. Running the pipeline with version 1.3.4 should provide the same outlier results as 1.3.3.
Due to snakemake updates affecting wBuild and the way we installed FRASER, installing DROP 1.3.3 no longer works. Version 1.3.4 simply fixes the FRASER version to ensure reproducibility and fixes certain scripts affected by the snakemake update. Running the pipeline with the version 1.3.4 should provide the same outlier results as 1.3.3.

Version 1.3.0 introduces the option to use FRASER 2.0 which is an improved version of FRASER that uses the Intron Jaccard Index metric instead of percent spliced in and splicing efficiency to quantify and later call aberrant splicing. To run FRASER 2.0, modify the `FRASER_version` parameter in the aberrantSplicing dictionary in the config file and adapt the `quantileForFiltering` and `deltaPsiCutoff` parameters. See the [config template](https://github.com/gagneurlab/drop/blob/master/drop/template/config.yaml) for more details. When switching between FRASER versions, we recommend running DROP in a
separate folder for each version. Moreover, DROP now allows users to provide lists of genes to focus on and do the multiple testing correction instead of the usual transcriptome-wide approach. Refer to the [documentation](https://gagneurlab-drop.readthedocs.io/en/latest/prepare.html#limiting-fdr-correction-to-subsets-of-genes-of-interest).

`Snakemake v.7.8` introduced some changes in which changes in parameters can cause rules to be re-executed. More info [here](https://github.com/snakemake/snakemake/issues/1694). This affects DROP and causes certain rules in the AS and QC modules to be triggered even if they were already completed and there were no changes in the sample annotation or scripts. The workaround is to run DROP by adding the parameter `--rerun-triggers mtime`, e.g. `snakemake -n --rerun-triggers mtime` or `snakemake --cores 10 --rerun-triggers mtime`. We will investigate the rules in DROP to fix this.

As of version 1.2.1 DROP has a new module that performs RNA-seq variant calling. The input are BAM files and the output either a single-sample or a multi-sample VCF file (option specified by the user) annotated with allele frequencies from gnomAD (if specified by the user). The sample annotation table does not need to be changed, but several new parameters in the config file have to be added and tuned. For more info, refer to the [documentation](https://gagneurlab-drop.readthedocs.io/en/latest/prepare.html#rna-variant-calling-dictionary).


Also, as of version 1.2.1 the integration of external split and non-split counts to detect aberrant splicing is now possible. In a new column in the sample annotation, simply specify the directory containing the counts. For more info, refer to the [documentation](https://gagneurlab-drop.readthedocs.io/en/latest/prepare.html#external-count-examples).



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
This shows you the rules of all subworkflows. Omit `-n` and specify the number of cores with `--cores ` if you are sure that you want you execute all printed rules. You can also invoke single workflows explicitly e.g. for aberrant expression with:
```
snakemake aberrantExpression --cores 10
```

## Datasets
The following publicly-available datasets of gene counts can be used as controls.
Please cite as instructed for each dataset.

* 154 non strand-specific fibroblasts, build hg19, Technical University of Munich: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4646822.svg)](https://doi.org/10.5281/zenodo.4646822)

* 135 strand-specific fibroblasts, build hg19, high seq depth (116 million mapped reads), Technical University of Munich: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7510836.svg)](https://zenodo.org/record/7510836)

* 127 strand-specific fibroblasts, build hg19, low seq depth (70 million mapped reads), Technical University of Munich: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7510845.svg)](https://zenodo.org/record/7510845)

* 49 tissues, each containing hundreds of samples, non strand-specific, build hg19, GTEx: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5596755.svg)](https://doi.org/10.5281/zenodo.5596755)

* 49 tissues, each containing hundreds of samples, non strand-specific, build hg38, GTEx: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6078396.svg)](https://doi.org/10.5281/zenodo.6078396)

* 139 strand-specific fibroblasts, build hg19, Baylor College of Medicine: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3963473.svg)](https://doi.org/10.5281/zenodo.3963473)

* 125 strand-specific blood, build hg19, Baylor College of Medicine: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3963470.svg)](https://doi.org/10.5281/zenodo.3963470)

* 330 strand-specific induced pluripotent stem cells (iPSCs), build hg19, EMBL: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7022459.svg)](https://doi.org/10.5281/zenodo.7022459)

* 56 non strand-specific amniotic fluid cells, build hg19, The University of Hong Kong: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7079684.svg)](https://doi.org/10.5281/zenodo.7079684)

If you want to contribute with your own count matrices, please contact us: yepez at in.tum.de

## Citation

If you use DROP in research, please cite our [manuscript](https://www.nature.com/articles/s41596-020-00462-5).

Furthermore, if you use:
* the aberrant expression module, also cite [OUTRIDER](https://doi.org/10.1016/j.ajhg.2018.10.025)
* the aberrant splicing module, also cite [FRASER](https://www.nature.com/articles/s41467-020-20573-7) or [FRASER2](https://www.sciencedirect.com/science/article/pii/S0002929723003671?dgcid=coauthor), depending on the version that you use
* the MAE module, also cite the [Kremer, Bader et al study](https://www.nature.com/articles/ncomms15824) and [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).

For the complete set of tools used by DROP (e.g. for counting), see the [manuscript](https://www.nature.com/articles/s41596-020-00462-5).

## Acknowledgements and Funding

The DROP team is composed of members from the Gagneur lab at the Department of Informatics and School of Medicine of the Technical University of Munich (TUM) and The German Human Genome-Phenome Archive (GHGA). The team has been funded by the German Bundesministerium für Bildung und Forschung (BMBF) through the e:Med Networking fonds AbCD-Net, Medical Informatics Initiative CORD-MI, and ERA PerMed project PerMiM. We would like to thank all the users for their feedback.
