Preparing the Input Data
========================

The input files of DROP are: 

- BAM files from RNA-seq  (and their respective index files)
- VCF files from either WES or WGS (and their respective index files). Only used for the MAE module
- a configuration file containing the different parameters
- a sample annotation file
- a gene annotation file (gtf)
- a reference genome file (fasta, and its respective index)

For more details see the Materials section of the `DROP manuscript <https://rdcu.be/cdMmF>`_.


Config file
-----------

The config file is in `YAML <https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html>`_ format. It is composed of general and module-specific parameters. In *YAML*, a variable can be of the following types: boolean, string, numeric, list and dictionary. They are declared by writing the variable name followed by a colon, a space, and the value, for example:

.. code-block:: yaml

    # A boolean is binary and can be true of false
    boolean_var: true    # or false

    # A string is a text. Quotation marks are not needed.
    string_var: whatever text

    # A numeric can be an integer or a real number. A dot separates a decimal.
    numeric_var: 0.05

    # A list is a collection of elements of the same type.
    list_var:
      - element_1   # elements are indented
      - element_2

    # A dictionary contains key-value pairs. It is a collection of multiple
    #    elements where the key is a string and the value any type.
    dictionary_var:
      key_1: value_1
      key_2: value_2


Now we describe the different parameters needed in DROP.
When providing a path to a file or directory, please provide the *full system path*. 

Global parameters
+++++++++++++++++
These parameters are applied to multiple modules and as a result should be consistent throughout the data you are analyzing

===================  ==========  =======================================================================================================================================  ======
Parameter            Type        Description                                                                                                                              Default/Examples
===================  ==========  =======================================================================================================================================  ======
projectTitle         character   Title of the project to be displayed on the rendered HTML output                                                                         ``Project 1``
htmlOutputPath       character   Full path of the folder where the HTML files are rendered                                                                                ``/data/project1/htmlOutput``
indexWithFolderName  boolean     If true, the basename of the project directory will be used as prefix for the index.html file                                            ``true``
genomeAssembly       character   Either hg19/hs37d5 or hg38/GRCh38, depending on the genome assembly used for mapping                                                     ``/data/project1``
sampleAnnotation     character   Full path of the sample annotation table                                                                                                 ``/data/project1/sample_annotation.tsv``
root                 character   Full path of the folder where the sub-directories processed_data and processed_results will be created containing DROP's output files.   ``/data/project1``
genome               character   Full path of a human reference genome fasta file                                                                                         ``/path/to/hg19.fa``
genome               dictionary  (Optional) Multiple fasta files can be specified when RNA-seq BAM files belong to different genome. assemblies (eg, ncbi, ucsc).         ``ncbi: /path/to/hg19_ncbi.fa``

                                                                                                                                                                          ``ucsc: /path/to/hg19_ucsc.fa``
geneAnnotation       dictionary  A key-value list of the annotation name (key) and the full path to the GTF file (value). More than one annotation file can be provided.  ``anno1: /path/to/gtf1.gtf``

                                                                                                                                                                          ``anno2: /path/to/gtf2.gtf``
hpoFile              character   Full path of the file containing HPO terms. If ``null`` (default), it reads it from our webserver. Refer to `files-to-download`_         ``/path/to/hpo_file.tsv``
tools                dictionary  A key-value list of different commands (key) and the command (value) to run them                                                         ``gatkCmd: gatk``

                                                                                                                                                                          ``bcftoolsCmd: bcftools``

                                                                                                                                                                          ``samtoolsCmd: samtools``
===================  ==========  =======================================================================================================================================  ======

Export counts dictionary
++++++++++++++++++++++++
These parameters are directly used by the ``exportCounts`` snakemake command. This section
is used to designate which aberrant expression and aberrant splicing groups should be exported
into datasets that can be shared. To avoid sharing sensitive data, only the canonical annotations
as described by `geneAnnotations` are exported. Only the groups excluded by `excludeGroups` are not exported.

===============  ====  ==========================================================================================================================  ======
Parameter        Type  Description                                                                                                                 Default/Examples
===============  ====  ==========================================================================================================================  ======
geneAnnotations  list  key(s) from the ``geneAnnotation`` parameter, whose counts should be exported                                               ``- gencode34``
excludeGroups    list  aberrant expression and aberrant splicing groups whose counts should not be exported. If ``null`` all groups are exported.  ``- group1``
===============  ====  ==========================================================================================================================  ======


Aberrant expression dictionary
++++++++++++++++++++++++++++++
These parameters are directly used by the ``aberrantExpression`` snakemake command. Aberrant expression groups must have at least ``10``
samples per group. To use external counts please see the ``Using External Counts`` section.

============================  =========  =======================================================================================================================================================================  ======
Parameter                     Type       Description                                                                                                                                                              Default/Examples
============================  =========  =======================================================================================================================================================================  ======
run                           boolean    If true, the module will be run. If false, it will be ignored.                                                                                                           ``true``
groups                        list       DROP groups that should be executed in this module. If not specified or ``null`` all groups are used.                                                                    ``- group1``

                                                                                                                                                                                                                  ``- group2``
minIds                        numeric    A positive number indicating the minimum number of samples that a group needs in order to be analyzed. We recommend at least 50.                                         ``1``
fpkmCutoff                    numeric    A positive number indicating the minimum FPKM per gene that 5% of the samples should have. If a gene has less it is filtered out.                                        ``1 # suggested by OUTRIDER``
implementation                character  Either 'autoencoder', 'pca' or 'peer'. Methods to remove sample covariation in OUTRIDER.                                                                                 ``autoencoder``
zScoreCutoff                  numeric    A non-negative number. Z scores (in absolute value) greater than this cutoff are considered as outliers.                                                                 ``0``
padjCutoff                    numeric    A number between (0, 1] indicating the maximum FDR an event can have in order to be considered an outlier.                                                               ``0.05``
maxTestedDimensionProportion  numeric    An integer that controls the maximum value that the encoding dimension can take. Refer to `advanced-options`_.                                                           ``3``
yieldSize                     numeric    An integer that sets the batch size for counting reads within a bam file. If memory issues persist lower the yieldSize.                                                  ``2000000``
genesToTest                   character  Full path to a yaml file specifying lists of candidate genes per sample to test during FDR correction. See the documentation for details on the structure of this file.  ``/path/to/genes_to_test.yaml``
============================  =========  =======================================================================================================================================================================  ======

Aberrant splicing dictionary
++++++++++++++++++++++++++++
These parameters are directly used by the ``aberrantSplicing`` snakemake command. Each group must have at least ``10``
samples. This module uses FRASER to detect aberrant splicing. We recently developed an improved version of FRASER that uses
the Intron Jaccard Index instead of percent spliced in and splicing efficiency to call aberrant splicing. To use this improved version,
set the ``FRASER_version`` parameter to 'FRASER2'. When switching between FRASER versions, we recommend running DROP in a
separate folder for each version.
To use external counts, refer to the ``Using External Counts`` section. 

============================  =========  =====================================================================================================================================================================================================================  ======
Parameter                     Type       Description                                                                                                                                                                                                            Default/Examples
============================  =========  =====================================================================================================================================================================================================================  ======
run                           boolean    If true, the module will be run. If false, it will be ignored.                                                                                                                                                         ``true``
groups                        list       Same as in aberrant expression.                                                                                                                                                                                        ``# see aberrant expression example``
minIds                        numeric    Same as in aberrant expression.                                                                                                                                                                                        ``1``
recount                       boolean    If true, it forces samples to be recounted.                                                                                                                                                                            ``false``
longRead                      boolean    Set to true only if counting Nanopore or PacBio long reads.                                                                                                                                                            ``false``
keepNonStandardChrs           boolean    Set to true if non standard chromosomes are to be kept for further analysis.                                                                                                                                           ``false``
filter                        boolean    If false, no filter is applied. We recommend filtering.                                                                                                                                                                ``true``
minExpressionInOneSample      numeric    The minimal read count in at least one sample required for an intron to pass the filter.                                                                                                                               ``20``
quantileMinExpression         numeric    The minimum total read count (N) an intron needs to have at the specified quantile across samples to pass the filter. See ``quantileForFiltering``.                                                                    ``10``
quantileForFiltering          numeric    Defines at which percentile the ``quantileMinExpression`` filter is applied. A value of 0.95 means that at least 5% of the samples need to have a total read count N >= ``quantileMinExpression`` to pass the filter.  ``0.95``
minDeltaPsi                   numeric    The minimal variation (in delta psi) required for an intron to pass the filter.                                                                                                                                        ``0.05``
implementation                character  Either 'PCA' or 'PCA-BB-Decoder'. Methods to remove sample covariation in FRASER.                                                                                                                                      ``PCA``
FRASER_version                character  Either 'FRASER' or 'FRASER2'.                                                                                                                                                                                          ``FRASER``
deltaPsiCutoff                numeric    A non-negative number. Delta psi values greater than this cutoff are considered as outliers. Set to 0.1 when using FRASER2.                                                                                            ``0.3 # suggested by FRASER``
padjCutoff                    numeric    Same as in aberrant expression.                                                                                                                                                                                        ``0.1``
maxTestedDimensionProportion  numeric    Same as in aberrant expression.                                                                                                                                                                                        ``6``
genesToTest                   character  Same as in aberrant expression.                                                                                                                                                                                        ``/path/to/genes_to_test.yaml``
============================  =========  =====================================================================================================================================================================================================================  ======


Mono-allelic expression (MAE) dictionary
++++++++++++++++++++++++++++++++++++++++
These parameters are directly used by the ``mae`` snakemake command. MAE groups are not bound by a minimum number of samples,
but require additional information in the sample annotation table.

=====================  =========  ========================================================================================================================  ======
Parameter              Type       Description                                                                                                               Default/Examples
=====================  =========  ========================================================================================================================  ======
run                    boolean    If true, the module will be run. If false, it will be ignored.                                                            ``true``
groups                 list       Same as in aberrant expression.                                                                                           ``# see aberrant expression example``
gatkIgnoreHeaderCheck  boolean    If true (recommended), it ignores the header warnings of a VCF file when performing the allelic counts                    ``true``
padjCutoff             numeric    Same as in aberrant expression.                                                                                           ``0.05``
allelicRatioCutoff     numeric    A number between [0.5, 1) indicating the maximum allelic ratio allele1/(allele1+allele2) for the test to be significant.  ``0.8``
addAF                  boolean    Whether or not to add the allele frequencies from gnomAD                                                                  ``true``
maxAF                  numeric    Maximum allele frequency (of the minor allele) cut-off. Variants with AF equal or below this number are considered rare.  ``0.001``
maxVarFreqCohort       numeric    Maximum variant frequency among the cohort.                                                                               ``0.05``
qcVcf                  character  Full path to the vcf file used for VCF-BAM matching. Refer to `files-to-download`_.                                       ``/path/to/qc_vcf.vcf.gz``
qcGroups               list       Same as “groups”, but for the VCF-BAM matching                                                                            ``# see aberrant expression example``
dnaRnaMatchCutoff      numeric    fraction (0-1) used to seperate "matching" samples and "non-matching" samples comparing the DNA and RNA data during QC    ``0.85``
=====================  =========  ========================================================================================================================  ======


RNA Variant Calling dictionary
++++++++++++++++++++++++++++++++++
Calling variants on RNA-seq data may be useful for researchers who do not have access to variant calls from genomic data. While variant calling from WES and WGS technologies may be more traditional (and reliable), variant calling from RNA-Seq data can provide additional evidence for the underlying causes of aberrant expression or splicing.
The RNA variant calling process uses information from multiple samples (as designated by the ``groups`` variable) to improve the quality of the called variants. However, the larger the group size, the more costly the computation is in terms of time and resources. To prioritize accuracy, include many samples in each ``DROP_GROUP``, and to prioritize speed up computation, separate samples into many groups. Additionally, certain vcf and bed files must be included to further boost the quality of the called variants (refer to `files-to-download`_).

=====================  =========  ================================================================================================================================================================================================  =========
Parameter              Type       Description                                                                                                                                                                                       Default/Examples
=====================  =========  ================================================================================================================================================================================================  =========
run                    boolean    If true, the module will be run. If false, it will be ignored.                                                                                                                                    ``true``
groups                 list       Same as in aberrant expression.                                                                                                                                                                   ``# see aberrant expression example``
highQualityVCFs        list       File paths where each item is the path to a vcf file. Each vcf file describes known high quality variants, which are used to recalibrate sequencing scores. Refer to `files-to-download`_         ``- known_indels.vcf``

                                                                                                                                                                                                                                    ``- known_SNPs.vcf``

dbSNP                  character  Location of the dbSNP ``.vcf`` file. This improves both recalibrating sequencing scores, as well as variant calling precision. Refer to `files-to-download`_                                      ``path/to/dbSNP.vcf``
repeat_mask            character  Location of the RepeatMask ``.bed`` file. Refer to `files-to-download`_                                                                                                                           ``path/to/RepeatMask.bed``
createSingleVCF        boolean    If ``true``, splits the multi-sample VCF file into individual sample VCF files. This only subsets the larger vcf sample.                                                                          ``true``
addAF                  boolean    Whether or not to add the allele frequencies from gnomAD                                                                                                                                          ``true``
maxAF                  numeric    Maximum allele frequency (of the minor allele) cut-off. Variants with AF equal or below this number are considered rare.                                                                          ``0.001``
maxVarFreqCohort       numeric    Maximum variant frequency among the cohort.                                                                                                                                                       ``0.05``
minAlt                 numeric    Integer describing the minimum required reads that support the alternative allele. We recommend a minimum of 3 if further filtering on your own. 10 otherwise.                                    ``3``
hcArgs                 character  String describing additional arguments for GATK haplocaller. Refer to `advanced-options`_.                                                                                                        ``""``
yieldSize              numeric    An integer that sets the batch size for counting reads within a vcf file. If memory issues persist during ``batch_data_table`` lower the yieldSize.                                               ``100000``
=====================  =========  ================================================================================================================================================================================================  =========


Modularization of DROP
-----------------------------------
DROP allows to control which modules to run via the  ``run`` variable in the config file. By default, each module is set to ``run: true``.  Setting this value to  ``false``  stops a particular module from being run. This will be noted as a warning at the beginning of the ``snakemake`` run, and the corresponding module will be renamed in the ``Scripts/`` directory.

For example, if the AberrantExpression module is set to false, the  ``Scripts/AberrantExpression/`` directory will be renamed to ``Scripts/_AberrantExpression/`` which tells DROP not to execute this module.


Creating the sample annotation table
------------------------------------
For a detailed explanation of the columns of the sample annotation, please refer to
Box 3 of the `DROP manuscript <https://rdcu.be/cdMmF>`_. Although some information has been updated since puplication, please use this documentation as the preferred syntax/formatting.

Each row of the sample annotation table corresponds to a unique pair of RNA and DNA
samples derived from the same individual. An RNA assay can belong to one or more DNA
assays, and vice-versa. If so, they must be specified in different rows. The required
columns are ``RNA_ID``, ``RNA_BAM_FILE`` and ``DROP_GROUP``, plus other module-specific
ones (see DROP manuscript).

The following columns describe the RNA-seq experimental setup:
``PAIRED_END``, ``STRAND``, ``COUNT_MODE`` and ``COUNT_OVERLAPS``. They affect the
counting procedures of the aberrant expression and splicing modules. For a detailed
explanation, refer to the documentation of `HTSeq <https://htseq.readthedocs.io/en/latest/>`_.

To run the MAE module, the columns ``DNA_ID`` and ``DNA_VCF_FILE`` are needed. MAE can not be run
in samples using external counts as we need to use the ``RNA_BAM_FILE`` to count reads supporting
each allele of the heterozygous variants found in the ``DNA_VCF_FILE``.

In case RNA-seq BAM files belong to different genome assemblies (eg, ncbi, ucsc), multiple
reference genome fasta files can be specified. Add a column called `GENOME` that
contains, for each sample, the key from the `genome` parameter in the config file that
matches its genome assembly (eg, ncbi or ucsc).

The sample annotation file must be saved in the tab-separated values (tsv) format. The 
column order does not matter. Also, it does not matter where it is stored, as the path is 
specified in the config file. Here we provide some examples on how to deal with certain
situations. For simplicity, we do not include all possible columns in the examples.

=====================  =========  ================================================================================================================================================================================================  ==========================
Parameter              Type       Description                                                                                                                                                                                       Default/Examples
=====================  =========  ================================================================================================================================================================================================  ==========================
RNA_ID                 character  Unique identifier from an RNA assay.                                                                                                                                                              ``sample1``
RNA_BAM_FILE           character  Absolute path of the BAM file derived from RNA-seq. A BAM file can belong to only one RNA_ID and vice versa.                                                                                      ``path/to/sample1.bam``                                                                          
DNA_VCF_FILE           character  Absolute path to the corresponding VCF. The DNA_ID has to match the ID inside the VCF file. In case a multisample VCF is used, write the file name for each sample.                               ``path/to/sample1.vcf``
DNA_ID                 character  Unique identifier from a DNA assay.                                                                                                                                                               ``sample1``
DROP_GROUP             list       The analysis group(s) that the RNA assay belongs to. Multiple groups must be separated by commas and no spaces (e.g. blood,WES,groupA). We recommend doing a different analysis for each tissue 
                                  as gene expression and splicing can be tissue specific.                                                                                                                                           ``group1,group2``
PAIRED_END             boolean    Either TRUE or FALSE, depending on whether the sample comes from paired-end RNA-seq or not.                                                                                                       ``TRUE``
COUNT_MODE             character  Either ``Union``, ``IntersectionStrict`` or ``IntersectionNotEmpty``. Refer to the documentation of HTSeq for details.                                                                            ``IntersectionStrict``
COUNT_OVERLAPS         character  Either TRUE or FALSE, depending on whether reads overlapping different regions are allowed and counted.                                                                                           ``TRUE``
STRAND                 character  Either yes, no, or reverse: ``no`` means that the sequencing was not strand specific; ``yes`` that it was strand specific, and the first read in the pair is on the same strand as the feature 
                                  and the second read on the opposite strand; and ``reverse`` that the sequencing is strand specific and the first read in the pair is on the opposite strand to the feature and the second read 
                                  on the same strand.                                                                                                                                                                               ``no``
HPO_TERMS              list       Comma-separated phenotypes encoded as HPO terms.                                                                                                                                                  ``HP:0001479, HP:0005591``
GENE_COUNTS_FILE       character  (Only required for aberrant expression external samples) Location of external gene-level count matrix.                                                                                            ``/path/to/gene_counts/``
GENE_ANNOTATION        character  (Only required for aberrant expression external samples) Gene annotation used to obtain the count matrix. Must correspond to the key of an entry in the geneAnnotation parameter of the config 
                       file.                                                                                                                                                                                                        ``v29``
GENOME                 character  (Optional) Either ``ncbi`` or ``ucsc`` indicating the reference genome assembly.                                                                                                                  ``ncbi``
SPLICE_COUNTS_DIR      character  (Only required for aberrant splicing external samples) Location of external files required for aberrant splicing module as explained above.                                                       ``/path/to/splicing_dir/``
SEX                    character  (Optional) Either ``m``, ``male``, ``f`` or ``female`` or ``unknown`` . When provided, sex matching algorithm will be run to match provided sex values to bam files and predict SEX value for 
                                  unknown samples.                                                                                                                                                                                  ``m``
TISSUE                 character  (Optional)                                                                                                                                                                                        ``BRAIN``
DISEASE                character  (Optional)                                                                                                                                                                                        ``AML``
=====================  =========  ================================================================================================================================================================================================  ==========================



Using External Counts
++++++++++++++++++++++++++++++++++
DROP can utilize external counts for the ``aberrantExpression`` and ``aberrantSplicing`` modules
which can enhance the statistical power of these modules by providing more samples from which we 
can build a distribution of counts and detect outliers. However this process introduces some
particular issues that need to be addressed to make sure it is a valuable addition to the experiment.

In case external counts are included, add a new row for each sample from those 
files (or a subset if not all samples are needed). Add the columns: ``GENE_COUNTS_FILE``
(for aberrant expression), ``GENE_ANNOTATON``, and ``SPLICE_COUNTS_DIR`` (for aberrant splicing).
These columns should remain empty for samples processed locally (from ``RNA_BAM``).

Aberrant Expression
####################
Using external counts for aberrant expression forces you to use the exact same gene annotation for each
external sample as well as using the same gene annotation file specified in the config file
``Global parameters`` section. This is to avoid potential mismatching on counting, 2 different gene
annotations could drastically affect which reads are counted in which region drastically skewing the results.

The user must also use special consideration when building the sample annotation table. Samples
using external counts need only ``RNA_ID`` which must exactly match the column header in the external count file
``DROP_GROUP``, ``GENE_COUNTS_FILE``, and ``GENE_ANNOTATION`` which must contain the exact key specified in the config.
The other columns should remain empty. 

Using ``exportCounts`` generates the sharable ``GENE_COUNTS_FILE`` file in the appropriate
``ROOT_DIR/Output/processed_results/exported_counts/`` sub-directory.

Aberrant Splicing
##################
Using external counts for aberrant splicing reduces the number of introns processed to only those
that are exactly the same between the local and external junctions. Because rare junctions may be 
personally identifiable the ``exportCounts`` command only exports regions canonically mentioned in the gtf file.
As a result, when merging the external counts with the local counts we only match introns that are **exact** between
the 2 sets, this is to ensure that if a region is missing we don't introduce 0 counts into the distribution calculations.

The user must also use special consideration when building the sample annotation table. Samples
using external counts need only ``RNA_ID`` which must exactly match the column header in the external count file
``DROP_GROUP``, and ``SPLICE_COUNTS_DIR``. ``SPLICE_COUNTS_DIR`` is the directory containing the set of 5 needed count files.
The other columns should remain empty. 

Using ``exportCounts`` generates the necessary files in the appropriate
``ROOT_DIR/Output/processed_results/exported_counts/`` sub-directory

``SPLICE_COUNTS_DIR`` should contain the following:  

* k_j_counts.tsv.gz  
* k_theta_counts.tsv.gz  
* n_psi3_counts.tsv.gz  
* n_psi5_counts.tsv.gz  
* n_theta_counts.tsv.gz  

Publicly available DROP external counts
#######################################
You can find different sets of publicly available external counts to add to your
analysis on our `github page <https://github.com/gagneurlab/drop/#datasets>`_

If you want to contribute with your own count matrices, please contact us: yepez at in.tum.de.

External count examples
+++++++++++++++++++++++

In case counts from external matrices are to be integrated into the analysis,
the sample annotation must be built in a particular way
A new row must be added for each sample from the count matrix that should be included in the 
analysis. The ``RNA_ID`` must match the column header of the external files,
the ``RNA_BAM_FILE`` must not be specified. The ``DROP_GROUP`` of the local
and external samples that are to be analyzed together must be the same.
For aberrant expression, the GENE_ANNOTATION of the external counts and the key of the `geneAnnotation`
parameter from the config file must match.

This example will use the ``DROP_GROUP`` BLOOD_AE for the aberrant expression module (containing S10R, EXT-1R, EXT-2R) and
the ``DROP_GROUP`` BLOOD_AS for the aberrant expression module (containing S10R, EXT-2R, EXT-3R)

======  ======  =================  =================  ==============================  =============== =========================
RNA_ID  DNA_ID  DROP_GROUP         RNA_BAM_FILE       GENE_COUNTS_FILE                GENE_ANNOTATION SPLICE_COUNTS_DIR
======  ======  =================  =================  ==============================  =============== =========================
S10R    S10G    BLOOD_AE,BLOOD_AS  /path/to/S10R.BAM  
EXT-1R          BLOOD_AE                              /path/to/externalCounts.tsv.gz  gencode34
EXT-2R          BLOOD_AE,BLOOD_AS                     /path/to/externalCounts.tsv.gz  gencode34       /path/to/externalCountDir 
EXT-3R          BLOOD_AS                                                                              /path/to/externalCountDir 
======  ======  =================  =================  ==============================  =============== =========================


Limiting FDR correction to subsets of genes of interest
-------------------------------------------------------
In addition to returning transcriptome-wide results, DROP provides the option to 
limit the FDR correction to user-provided genes of interest in the 
``aberrantExpression`` and ``aberrantSplicing`` modules. These could, for example, be all 
OMIM genes. It is also possible to provide sample-specific genes such as all 
genes with a rare splice-region variant for each sample. 
To use this feature, a YAML file containing the set(s) of genes to test 
(per sample or for all samples) needs to be specified in the ``genesToTest`` parameter 
of the ``aberrantExpression`` and ``aberrantSplicing`` modules in the config file. 
If no file is provided, only transcriptome-wide results will be reported.
Otherwise, the results tables of the ``aberrantExpression`` and ``aberrantSplicing`` modules 
will additionally report aberrant events passing the cutoffs based on calculating 
the FDR with respect to the genes in the provided lists.

Creating the YAML file specifying subsets of genes to test
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The file containing the list of genes (HGNC symbols) to be tested must be a YAML file, 
where the variable names specify the name of each set of tested genes. In the output 
of DROP, this name will be used to identify the set in the results table. Each set 
can either be: i) a list of genes, in which case the set will be tested for all samples, or ii)
sample-specific sets that can be created by giving the RNA_ID of the sample
for which the set should be used as the name (see example below).
This YAML file can be created in R using ``yaml::write_yaml(subsetList, filepath)``, 
where ``subsetList`` is a named list of named lists containing the sets of genes to test.
The gene names must match those from the provided gtf file. We currently do not support Ensembl ids as input.
The table with extracted gene names from the gtf file is located under: 
``root/processed_data/preprocess/{geneAnnotation}/gene_name_mapping_{geneAnnotation}.tsv``.
In the following example, the name of the global set of genes is ``Genes_to_test_on_all_samples``
and the name of the sample-specific set is ``Genes_with_rare_splice_variants``:

Example content of ``/path/to/genes_to_test.yaml``:

.. code-block:: bash

    Genes_to_test_on_all_samples:
      - BTG3
      - GATD3B
      - PKNOX1
      - APP
      - RRP1
      - WRB-SH3BGR
      - SLC19A1
    Genes_with_rare_splice_variants:
      sample1:
      - ABCG1
      - MCOLN1
      - SLC45A1
      sample2:
      - CLIC6
      - ATP5PO
      - WRB
      - ETS2
      - HLCS



.. _files-to-download:

Files to download
-----------------

The following files can be downloaded from our `public repository <https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/>`_.

1. VCF file containing different positions to be used to match DNA with RNA files.
The file name is ``qc_vcf_1000G_{genome_build}.vcf.gz``. One file is available for each
genome build (hg19/hs37d5 and hg38/GRCh38). Download it together with the corresponding ``.tbi`` file.
Indicate the full path to the vcf file in the ``qcVcf`` key in the mono-allelic expression dictionary.
This file is only needed for the MAE module. Otherwise, write ``null`` in the ``qcVcf`` key.

2. Text file containing the relations between genes and phenotypes encoded as HPO terms.
The file name is ``hpo_genes.tsv.gz``.
Download it and indicate the full path to it in the ``hpoFile`` key.
The file is only needed in case HPO terms are specified in the sample annotation.
Otherwise, write ``null`` in the ``hpoFile`` key.

3. For the ``rnaVariantCalling`` module known high quality variants are needed to calibrate variant and sequencing scores to be used in the ``rnaVariantCalling`` module in the ``highQualityVCF`` config parameter.
These and the associated ``.tbi`` indexes can be downloaded for hg19 at our `public repository <https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/>`_
and for hg38 through the Broad Institute's `resource bundle. <https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle>`_

hg19

* ``Mills_and_1000G_gold_standard.indels.hg19.sites.chrPrefix.vcf.gz``
* ``1000G_phase1.snps.high_confidence.hg19.sites.chrPrefix.vcf.gz``

hg38

* ``Mills_and_1000G_gold_standard.indels.hg38.vcf.gz``
* ``Homo_sapiens_assembly38.known_indels.vcf.gz``

We also recommend using the variants from dbSNP which is quite large. You can download them and their associated ``.tbi`` indexes from `NCBI <https://ftp.ncbi.nih.gov/snp/organisms/>`_

* follow links for the current version (``human_9606/VCF/00-All.vcf.gz``) or older assemblies (eg. ``human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz``)

The repeat masker file is used to filter hard-to-call regions. In general, this removes false-positive calls, however, some targeted and known splicing defects lie within these repeat regions. Understand that this filter is labeled ``Mask`` in the result VCF files. You can download the repeat mask and associated ``.idx`` on our `public repository. <https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/>`_ for the ``repeat_mask`` config parameter.

Example of RNA replicates 
-------------------------

======  ======  ==========  ===================  ==
RNA_ID  DNA_ID  DROP_GROUP  RNA_BAM_FILE         DNA_VCF_FILE
======  ======  ==========  ===================  ==
S10R_B  S10G    BLOOD       /path/to/S10R_B.BAM  /path/to/S10G.vcf.gz
S10R_M  S10G    MUSCLE      /path/to/S10R_M.BAM  /path/to/S10G.vcf.gz
======  ======  ==========  ===================  ==

Example of DNA replicates
++++++++++++++++++++++++++++++++++

======  ======  ==========  =================  ==
RNA_ID  DNA_ID  DROP_GROUP  RNA_BAM_FILE       DNA_VCF_FILE
======  ======  ==========  =================  ==
S20R    S20E    WES         /path/to/S20R.BAM  /path/to/S20E.vcf.gz
S20R    S20G    WGS         /path/to/S20R.BAM  /path/to/S20G.vcf.gz
======  ======  ==========  =================  ==

Example of a multi-sample vcf file
++++++++++++++++++++++++++++++++++

======  ======  ==========  =================  ==
RNA_ID  DNA_ID  DROP_GROUP  RNA_BAM_FILE       DNA_VCF_FILE
======  ======  ==========  =================  ==
S10R    S10G    WGS         /path/to/S10R.BAM  /path/to/multi_sample.vcf.gz
S20R    S20G    WGS         /path/to/S20R.BAM  /path/to/multi_sample.vcf.gz
======  ======  ==========  =================  ==

.. _advanced-options:

Advanced options
----------------

A local copy of DROP can be edited and modified.
For example, the user might want to add new plots to the ``Summary`` scripts, add
additional columns to the results tables, or modify the number of threads allowed for a script.

.. note::

    DROP needs to be installed from a local directory using ``pip install -e <path/to/drop-repo>``
    so that any changes in the code will be available in the next pipeline run
    Any changes made to the R code need to be updated with ``drop update`` in the project directory.

The aberrant expression and splicing modules use a denoising autoencoder to
correct for sample covariation. This process reduces the fitting space to a
dimension smaller than the number of samples N. The encoding dimension is optimized.
By default, the maximum value in the search space is N/3 for the aberrant expression,
and N/6 for the aberrant splicing case. The user can specify the
denominator with the parameter ``maxTestedDimensionProportion``.

DROP allows that BAM files from RNA-seq from samples belonging to the same `DROP_GROUP`
were aligned to different genome assemblies from the same build (e.g., some to ucsc
and others to ncbi, but all to either hg19 or hg38). If so, for the aberrant
expression and splicing modules, no special configuration is needed.
For the MAE and rnaVariantCalling modules, the different fasta files must be specified as a dictionary in
the `genome` parameter of the config file, and, for each sample, the corresponding
key of the `genome` dictionary must be specified in the `GENOME` column of the
sample annotation.
In additon, DROP allows that BAM files from RNA-seq were aligned to one genome
assembly (eg ucsc) and the corresponding VCF files from DNA sequencing to another
genome assembly (eg ncbi). If so, the assembly of the reference genome fasta file
must correspond to the one of the BAM file from RNA-seq.

Specific haplotype parameters can be denoted in the config file to further customize the RNA-seq variant calling. 
The different available parameters can be found in the
`HaplotypeCaller GATK documentation. <https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller>`_
One example for the value in the config file would be "--assembly-region-padding 100 --base-quality-score-threshold 18".
