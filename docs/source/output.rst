Results and Output of DROP
===========================

DROP is intended to help researchers use RNA-Seq data in order to detect genes with aberrant expression,
aberrant splicing and mono-allelic expression. By simplifying the workflow process we hope to provide
easy-to-read HTML files and output files. This section explains the results files. The paths of the output
files correspond to the ones from the demo (that can be run with the following code snippet)::

    #install drop
    mamba create -n drop_env -c conda-forge -c bioconda drop
    conda activate drop_env
    
    mkdir drop_demo
    cd drop_demo
    drop demo
    
    snakemake -c1

Aberrant Expression
+++++++++++++++++++

HTML file
#########
Looking at the resulting ``Output/html/drop_demo_index.html`` we can see the ``AberrantExpression`` 
tab at the top of the screen. The Overview tab contains links to the:  

* Counts Summaries for each aberrant expression group
    * number of local and external samples
    * Mapped reads and size factors for each sample
    * histograms showing the mean count distribution with different conditions
    * expressed genes within each sample and as a dataset
* Outrider Summaries for each aberrant expression group
    * aberrantly expressed genes per sample
    * correlation between samples before and after the autoencoder
    * biological coefficient of variation
    * aberrant samples
    * results table
* Files for each aberrant expression group
    * OUTRIDER datasets 
        * Follow the `OUTRIDER vignette <https://www.bioconductor.org/packages/devel/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf>`_ for individual OUTRIDER object file (ods) analysis.
    * Results tables
        * ``results.tsv`` this text file contains only the significant genes and samples that meet the cutoffs defined in the config file for ``padjCutoff`` and ``zScoreCutoff``

Local result files
##################
Additionally the ``aberrantExpression`` module creates the file ``Output/processed_results/aberrant_expression/{annotation}/outrider/{drop_group}/OUTRIDER_results_all.Rds``. This file contains the entire OUTRIDER results table regardless of significance.

Aberrant Splicing
+++++++++++++++++

HTML file
##########
Looking at the resulting ``Output/html/drop_demo_index.html`` we can see the ``AberrantSplicing`` 
tab at the top of the screen. The Overview tab contains links to the:  

* Counting Summaries for each aberrant splicing group
    * number of local and external samples
    * number introns/splice sites before and after merging
    * comparison of local and external mean counts
    * histograms showing the junction expression before and after filtering and variability
* FRASER Summaries for each aberrant splicing group
    * the number of samples, introns, and splice sites 
    * correlation between samples before and after the autoencoder
    * results table
* Files for each aberrant splicing group
    * FRASER datasets (fds)
        * Follow the `FRASER vignette <https://www.bioconductor.org/packages/devel/bioc/vignettes/FRASER/inst/doc/FRASER.pdf>`_ for individual FRASER object file (fds) analysis.
    * Results tables
        * ``results_per_junction.tsv`` this text file contains only significant junctions that meet the cutoffs defined in the config file. 

Local result files
##################
Additionally the ``aberrantSplicing`` module creates the following file ``Output/processed_results/aberrant_splicing/results/{annotation}/fraser/{drop_group}/results.tsv``.
This text file contains only significant junctions that meet the cutoffs defined in the config file, aggregated at the gene level. Any sample/gene pair is represented by only the most significant junction.

Mono-allelic Expression
+++++++++++++++++++++++

HTML file
##########
Looking at the resulting ``Output/html/drop_demo_index.html`` we can see the ``MonoallelicExpression`` 
tab at the top of the screen. The Overview tab contains links to the:  

* Results for each mae group
    * number of samples, genes, and mono-allelically expressed heterozygous SNVs
    * a cascade plot that shows additional filters
    * histogram of inner cohort frequency
    * summary of the cascade plot and results table
* Files for each mae group
    * Allelic counts
        * a directory containing the allelic counts of heterozygous variants
    * Results data tables of each sample (.Rds)
        * Rds objects containing the full results table regardless of MAE status
    * Significant MAE results tables
        * a link to the results file
        * Only contains significant MAE for the alternative allele results and results that pass the config file cutoffs
* Quality Control
    * QC Overview
        * For each mae group QC checks for DNA/RNA matching
    
Local result files
##################
Additionally the ``mae`` module creates the following files:

* ``Output/processed_results/mae/{drop_group}/MAE_results_all_{annotation}.tsv.gz``
    * this file contains the MAE results of all heterozygous SNVs regardless of significance
* ``Output/processed_results/mae/{drop_group}/MAE_results_{annotation}.tsv``
    * this is the file linked in the HTML document and described above
* ``Output/processed_results/mae/{drop_group}/MAE_results_{annotation}_rare.tsv``
    * this file is a subset of ``MAE_results_{annotation}.tsv`` with only the variants that pass the allele frequency cutoffs. If ``add_AF`` is set to ``true`` in config file must meet minimum AF set by ``max_AF``. Additionally, the inner-cohort frequency must meet the ``maxVarFreqCohort`` cutoff

RNA Variant Calling
+++++++++++++++++++++++

HTML file
##########
Looking at the resulting ``Output/html/drop_demo_index.html`` we can see the ``rnaVariantCalling`` 
tab at the top of the screen. The Overview tab contains links to the:  

* Results for each rvc batch
    * a table summarizing the variants and genotypes that pass the variant calling filters for each sample
        * ``FILTER``: explained below
        * ``cohortFreq``: frequency of the variant within the batch (number of samples with the variant / total samples)
        * ``MAX_AF``: frequency of the variant from **gnomAD** if enabled
    * a subset table showing only the ``PASS_rare`` variants
* Boxplot and underyling table showing the distribution of variants and the effect of various filters, split by genotype
    * The following labels are applied and are not excluded from the final ``.vcf`` files.
        * ``PASS_common``: passes variant calling thresholds and fails either ``max_AF`` or ``maxVarFreqCohort`` cutoffs
        * ``PASS_rare``: passes variant calling thresholds and config ``max_AF`` and ``maxVarFreqCohort`` cutoffs
        * ``Seq_filter``: fails one of the default variant calling filters
        * ``Mask``: variant falls in a repeat/mask region
        * ``minALT``: variant passes ``Seq_filter`` but doesn't meet config ``minALT`` criteria
* Boxplot and underyling table showing the number of variants that pass or fail the different filters
    
Local result files
##################
Additionally the ``rnaVariantCalling`` module creates the following output directories:

* ``Output/processed_results/rnaVariantCalling/batch_vcfs``
    * this directory contains the multi-sample vcf files for each batch
* ``Output/processed_results/rnaVariantCalling/sample_vcfs``
    * this directory contains the single-sample vcf files if the config parameter ``createSingleVCF`` is set to ``true``
* ``Output/processed_results/rnaVariantCalling/data_tables``
    * this directory contains data tables for each batch of vcfs
