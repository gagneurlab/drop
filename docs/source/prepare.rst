Preparing the Input Data
========================

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

Global parameters
+++++++++++++++++

===================  ==========  =======================================================================================================================================  ======
Parameter            Type        Description                                                                                                                              Default/Examples
===================  ==========  =======================================================================================================================================  ======
projectTitle         character   Title of the project to be displayed on the rendered HTML output                                                                         ``Project 1``
htmlOutputPath       character   Full path of the folder where the HTML files are rendered                                                                                ``/data/project1/htmlOutput``
indexWithFolderName  boolean     If true, the basename of the project directory will be used as prefix for the index.html file                                            ``true``
genomeAssembly       character   Either hg19 or hg38, depending on the genome assembly used for mapping                                                                   ``/data/project1``
sampleAnnotation     character   Full path of the sample annotation table                                                                                                 ``/data/project1/sample_annotation.tsv``
root                 character   Full path of the folder where the subdirectories processed_data and processed_results will be created containing DROP's output files.    ``/data/project1``
geneAnnotation       dictionary  A key-value list of the annotation name (key) and the full path to the GTF file (value). More than one annotation file can be provided.  ``anno1: /path/to/gtf1.gtf``

                                                                                                                                                                          ``anno2: /path/to/gtf2.gtf``
hpoFile              character   Full path of the file containing HPO terms. If ``null`` (default), it reads it from our webserver. Refer to :ref:`filesdownload`.        ``/path/to/hpo_file.tsv``                                           
tools                dictionary  A key-value list of different commands (key) and the command (value) to run them                                                         ``gatkCmd: gatk``

                                                                                                                                                                          ``bcftoolsCmd: bcftools``

                                                                                                                                                                          ``samtoolsCmd: samtools``
===================  ==========  =======================================================================================================================================  ======

Export counts dictionary
++++++++++++++++++++++++

===============  ====  ==========================================================================================================================  ======
Parameter        Type  Description                                                                                                                 Default/Examples
===============  ====  ==========================================================================================================================  ======
geneAnnotations  list  key(s) from the ``geneAnnotation`` parameter, whose counts should be exported                                               ``- gencode34``
excludeGroups    list  aberrant expression and aberrant splicing groups whose counts should not be exported. If ``null`` all groups are exported.  ``- group1``
===============  ====  ==========================================================================================================================  ======


Aberrant expression dictionary
++++++++++++++++++++++++++++++

============================  =========  =================================================================================================================================  ======
Parameter                     Type       Description                                                                                                                        Default/Examples
============================  =========  =================================================================================================================================  ======
groups                        list       DROP groups that should be executed in this module. If not specified or ``null`` all groups are used.                              ``- group1``

                                                                                                                                                                            ``- group2``
minIds                        numeric    A positive number indicating the minimum number of samples that a group needs in order to be analyzed. We recommend at least 50.   ``1``
fpkmCutoff                    numeric    A positive number indicating the minimum FPKM 5% of the samples per gene should have. If a gene has less it will be filtered out.  ``1 # suggested by OUTRIDER``
implementation                character  Either 'autoencoder', 'pca' or 'peer'. Methods to remove sample covariation in OUTRIDER.                                           ``autoencoder``
zScoreCutoff                  numeric    A non-negative number. Z scores (in absolute value) greater than this cutoff are considered as outliers.                           ``0``
padjCutoff                    numeric    A number between (0, 1] indicating the maximum FDR an event can have in order to be considered an outlier.                         ``0.05``
maxTestedDimensionProportion  numeric    An integer that controls the maximum value that the encoding dimension can take. Refer to :ref:`advancedoptions`.                  ``3``
============================  =========  =================================================================================================================================  ======

Aberrant splicing dictionary
++++++++++++++++++++++++++++

============================  =========  ============================================================================================  ======
Parameter                     Type       Description                                                                                   Default/Examples
============================  =========  ============================================================================================  ======
groups                        list       Same as in aberrant expression.                                                               ``# see aberrant expression example``
minIds                        numeric    Same as in aberrant expression.                                                               ``1``
recount                       boolean    If true, it forces samples to be recounted.                                                   ``false``
longRead                      boolean    Set to true only if counting Nanopore or PacBio long reads.                                   ``false``
keepNonStandardChrs           boolean    Set to true if non standard chromosomes are to be kept for further analysis.                  ``true``                        
filter                        boolean    If false, no filter is applied. We recommend filtering.                                       ``true``
minExpressionInOneSample      numeric    The minimal read count in at least one sample required for an intron to pass the filter.      ``20``
minDeltaPsi                   numeric    The minimal variation (in delta psi) required for an intron to pass the filter.               ``0.05``
implementation                character  Either 'PCA' or 'PCA-BB-Decoder'. Methods to remove sample covariation in FRASER.             ``PCA``
deltaPsiCutoff                numeric    A non-negative number. Delta psi values greater than this cutoff are considered as outliers.  ``0.3 # suggested by FRASER``
padjCutoff                    numeric    Same as in aberrant expression.                                                               ``0.1``
maxTestedDimensionProportion  numeric    Same as in aberrant expression.                                                               ``6``
============================  =========  ============================================================================================  ======


Mono-allelic expression dictionary
++++++++++++++++++++++++++++++++++

=====================  =========  ========================================================================================================================  ======
Parameter              Type       Description                                                                                                               Default/Examples
=====================  =========  ========================================================================================================================  ======
groups                 list       Same as in aberrant expression.                                                                                           ``# see aberrant expression example``
genome                 character  Full path of a human reference genome fasta file                                                                          ``/path/to/hg19.fa``
gatkIgnoreHeaderCheck  boolean    If true (recommended), it ignores the header warnings of a VCF file when performing the allelic counts                    ``true``
padjCutoff             numeric    Same as in aberrant expression.                                                                                           ``0.05``
allelicRatioCutoff     numeric    A number between [0.5, 1) indicating the maximum allelic ratio allele1/(allele1+allele2) for the test to be significant.  ``0.8``
addAF                  boolean    Whether or not to add the allele frequencies from gnomAD                                                                  ``true``
maxAF                  numeric    Maximum allele frequency (of the minor allele) cut-off. Variants with AF equal or below this number are considered rare.  ``0.001``
maxVarFreqCohort       numeric    Maximum variant frequency among the cohort.                                                                               ``0.05``      
qcVcf                  character  Full path to the vcf file used for VCF-BAM matching. Refer to :ref:`filesdownload`.                                       ``/path/to/qc_vcf.vcf.gz``
qcGroups               list       Same as “groups”, but for the VCF-BAM matching                                                                            ``# see aberrant expression example``
=====================  =========  ========================================================================================================================  ======


Creating the sample annotation table
------------------------------------

For a detailed explanation of the columns of the sample annotation, please refer to
the DROP manuscript. 
Inside the sample annotation, each row corresponds to a unique pair of RNA and DNA
samples derived from the same individual. An RNA assay can belong to one or more DNA
assays, and vice-versa. If so, they must be specified in different rows. The required
columns are ``RNA_ID``, ``RNA_BAM_FILE`` and ``DROP_GROUP``, plus other module-specific
ones (see DROP manuscript). In case external counts are included, add a new row for each
sample from those files (or a subset if not all samples are needed).

The sample annotation file should be saved in the tab-separated values (tsv) format. The 
column order does not matter. Also, it does not matter where it is stored, as the path is 
specified in the config file. Here we provide some examples on how to deal with certain
situations. For simplicity, we do not include all possible columns in the examples.

Example of RNA replicates 
++++++++++++++++++++++++++++++++++

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

External count matrices
+++++++++++++++++++++++

In case counts from external matrices are to be integrated into the analysis,
the file must be specified in the GENE_COUNTS_FILE column. A new row must be
added for each sample from the count matrix that should be included in the 
analysis. An RNA_BAM_FILE must not be specified. The DROP_GROUP of the local
and external samples that are to be analyzed together must be the same.
Similarly, the GENE_ANNOTATION of the external counts and the key of the `geneAnnotation`
parameter from the config file must match.

======  ======  ==========  =================  ==============================  ==
RNA_ID  DNA_ID  DROP_GROUP  RNA_BAM_FILE       GENE_COUNTS_FILE                GENE_ANNOTATION
======  ======  ==========  =================  ==============================  ==
S10R    S10G    BLOOD       /path/to/S10R.BAM  
EXT-1R          BLOOD                          /path/to/externalCounts.tsv.gz  gencode34
EXT-2R          BLOOD                          /path/to/externalCounts.tsv.gz  gencode34
======  ======  ==========  =================  ==============================  ==

.. _filesdownload:

Files to download
-----------------

Two different files can be downloaded from our `public repository <https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/>`_. 

1) VCF file containing different positions to be used to match DNA with RNA files.
The file name is ``qc_vcf_1000G_{genome_build}.vcf.gz``. One file is available for each 
genome build (hg19 and hg38). Download it together with the corresponding .tbi file. 
Indicate the full path to the vcf file in the ``qcVcf`` key in the mono-allelic expression dictionary.
This file is only needed for the MAE module. Otherwise, write ``null`` in the 
``qcVcf``key.

2) Text file containing the relations between genes and phenotypes encoded as HPO terms. 
The file name is ``hpo_genes.tsv.gz``
Download it and indicate the full path to it in the ``hpoFile`` key.
The file is only needed in case HPO terms are specified in the sample annotation.
Otherwise, write ``null`` in the ``hpoFile`` key.

.. _advancedoptions:

Advanced options
----------------

A local copy of DROP can be edited and modified for uncovering potential issues or increasing outputs.
For example, the user might want to add new plots to the ``Summary`` scripts, or add
additional columns to the results tables.
Specifically, the number of threads allowed for a computational step can be modified by the user.

.. note::

    DROP needs to be installed from a local directory :ref:`otherversions` using ``pip install -e <path/to/drop-repo>``
    so that any changes in the code will be available in the next pipeline run.
    Any changes made to the R code need to be updated with ``drop update`` in the project directory.

The aberrant expression and splicing modules use a denoising autoencoder to
correct for sample covariation. This process reduces the fitting space to a 
dimension smaller than the number of samples N. The encoding dimension is optimized.
We recommend the search space to be at most N/3 for the aberrant expression, 
and N/6 for the aberrant splicing case. Nevertheless, the user can specify the 
denominator with the parameter ``maxTestedDimensionProportion``.


