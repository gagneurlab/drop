.. _prepare:

Preparing the Input Data
========================

Filling the Config File
-----------------------

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
indexWithFolderName  boolean     variable needed for wBuild, do not edit it                                                                                               ``true``
fileRegex            character   variable needed for wBuild, do not edit it                                                                                               ``.*\.R``
root                 character   Full path of the folder where the subdirectories processed_data and processed_results will be created containing DROP's output files.    ``/data/project1``
sampleAnnotation     character   Full path of the sample annotation table                                                                                                 ``/data/project1/sample_annotation.tsv``
geneAnnotation       dictionary  A key-value list of the annotation name (key) and the full path to the GTF file (value). More than one annotation file can be provided.  ``anno1: /path/to/gtf1.gtf``

                                                                                                                                                                          ``anno2: /path/to/gtf2.gtf``
tools                dictionary  A key-value list of different commands (key) and the command (value) to run them                                                         ``gatkCmd: gatk``

                                                                                                                                                                          ``bcftoolsCmd: bcftools``

                                                                                                                                                                          ``samtoolsCmd: samtools``
===================  ==========  =======================================================================================================================================  ======


Aberrant expression dictionary
++++++++++++++++++++++++++++++

================  =========  =====================================================================================================================================  ======
Parameter         Type       Description                                                                                                                            Default/Examples
================  =========  =====================================================================================================================================  ======
groups            list       DROP groups that should be executed in this module. If not specified, all groups are used.                                             ``- group1``

                                                                                                                                                                    ``- group2``
fpkmCutoff        numeric    A non-negative number indicating the minimum FPKM 5% of the samples per gene should have. If a gene has less it will be filtered out.  ``1 # suggested by OUTRIDER``
implementation    character  Either 'autoencoder', 'pca' or 'peer'. Methods to remove sample covariation in OUTRIDER.                                               ``autoencoder``
zScoreCutoff      numeric    A non-negative number. Z scores (in absolute value) greater than this cutoff are considered as outliers.                               ``0``
padjCutoff        numeric    A number between (0, 1] indicating the maximum FDR an event can have in order to be considered an outlier.                             ``0.05``
================  =========  =====================================================================================================================================  ======

Aberrant splicing dictionary
++++++++++++++++++++++++++++

========================  =========  ============================================================================================  ======
Parameter                 Type       Description                                                                                   Default/Examples
========================  =========  ============================================================================================  ======
groups                    list       Same as in aberrant expression.                                                               ``# see aberrant expression example``
recount                   boolean    If true, it forces samples to be recounted                                                    ``false``
longRead                  boolean    Set to true only if counting Nanopore or PacBio long reads.                                   ``false``
minExpressionInOneSample  numeric    The minimal read count in at least one sample required for an intron to pass the filter.      ``20``
correction                character  Either 'PCA' or 'PCA-BB-Decoder'. Methods to remove sample covariation in FRASER.             ``PCA``
deltaPsiCutoff            numeric    A non-negative number. Delta psi values greater than this cutoff are considered as outliers.  ``0.3 # suggested by FRASER``
padjCutoff                numeric    Same as in aberrant expression.                                                               ``0.1``
========================  =========  ============================================================================================  ======


Mono-allelic expression dictionary
++++++++++++++++++++++++++++++++++

==================  =========  ========================================================================================================================  ======
Parameter           Type       Description                                                                                                               Default/Examples
==================  =========  ========================================================================================================================  ======
groups              list       Same as in aberrant expression.                                                                                           ``# see aberrant expression example``
geneAssembly        character  Either hg19 or hg38, depending on the genome build used                                                                   ``hg19``
genome              character  Full path of a human reference genome fasta file                                                                          ``/path/to/hg19.fa``
padjCutoff          numeric    Same as in aberrant expression.                                                                                           ``0.05``
allelicRatioCutoff  numeric    A number between [0.5, 1) indicating the maximum allelic ratio allele1/(allele1+allele2) for the test to be significant.  ``0.8``
addAF               boolean    Whether or not to add the allele frequencies from gnomAD                                                                  ``true``
maxAF               numeric    Maximum allele frequency (of the minor allele) cut-off. Variants with AF equal or below this number are considered rare.  ``0.001``
qcVcf               character  Full path to the vcf file used for VCF-BAM matching                                                                       ``/path/to/qc_vcf.vcf.gz``
qcGroups            list       Same as “groups”, but for the VCF-BAM matching                                                                            ``# see aberrant expression example``
==================  =========  ========================================================================================================================  ======


Creating the Sample Annotation Table
------------------------------------

For details on how to generate the sample annotation, please refer to the DROP paper. Here we provide some examples.

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

======  ======  ==========  ===================  ==
RNA_ID  DNA_ID  DROP_GROUP  RNA_BAM_FILE         DNA_VCF_FILE
======  ======  ==========  ===================  ==
S20R    S20E    WES         /path/to/S20R.BAM    /path/to/S20E.vcf.gz
S20R    S20G    WGS         /path/to/S20R.BAM    /path/to/S20G.vcf.gz
======  ======  ==========  ===================  ==

Example of a multi-sample vcf file
++++++++++++++++++++++++++++++++++

======  ======  ==========  ===================  ==
RNA_ID  DNA_ID  DROP_GROUP  RNA_BAM_FILE         DNA_VCF_FILE
======  ======  ==========  ===================  ==
S10R    S10G    WGS         /path/to/S10R.BAM    /path/to/multi_sample.vcf.gz
S20R    S20G    WGS         /path/to/S20R.BAM    /path/to/multi_sample.vcf.gz
======  ======  ==========  ===================  ==
