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
genomeAssembly       character   Either hg19 or hg38, depending on the genome assembly used for mapping                                                                   ``/data/project1``
sampleAnnotation     character   Full path of the sample annotation table                                                                                                 ``/data/project1/sample_annotation.tsv``
root                 character   Full path of the folder where the subdirectories processed_data and processed_results will be created containing DROP's output files.    ``/data/project1``
geneAnnotation       dictionary  A key-value list of the annotation name (key) and the full path to the GTF file (value). More than one annotation file can be provided.  ``anno1: /path/to/gtf1.gtf``

                                                                                                                                                                          ``anno2: /path/to/gtf2.gtf``
scanBamParam         character   Either null or the path to an Rds file containing a scanBamParam object. Refer to the advanced options below.                            ``/path/to/scanBamParam.Rds``
tools                dictionary  A key-value list of different commands (key) and the command (value) to run them                                                         ``gatkCmd: gatk``

                                                                                                                                                                          ``bcftoolsCmd: bcftools``

                                                                                                                                                                          ``samtoolsCmd: samtools``
===================  ==========  =======================================================================================================================================  ======

exportCounts dictionary
++++++++++++++++++++++++++++++

===============  ====  ==========================================================================================================================  ======
Parameter        Type  Description                                                                                                                 Default/Examples
===============  ====  ==========================================================================================================================  ======
geneAnnotations  list  key(s) from the ``geneAnnotation`` parameter, whose counts should be exported                                               ``- v34``
excludeGroups    list  aberrant expression and aberrant splicing groups whose counts should not be exported. If ``null`` all groups are exported.  ``- group1``
===============  ====  ==========================================================================================================================  ======


Aberrant expression dictionary
++++++++++++++++++++++++++++++

============================  =========  =====================================================================================================================================  ======
Parameter                     Type       Description                                                                                                                            Default/Examples
============================  =========  =====================================================================================================================================  ======
groups                        list       DROP groups that should be executed in this module. If not specified or ``null`` all groups are used.                                  ``- group1``

                                                                                                                                                                                ``- group2``
minIds                        numeric    A non-negative number indicating the minimum number of samples that a group needs in order to be analyzed. We recommend at least 50.   ``1``
fpkmCutoff                    numeric    A non-negative number indicating the minimum FPKM 5% of the samples per gene should have. If a gene has less it will be filtered out.  ``1 # suggested by OUTRIDER``
implementation                character  Either 'autoencoder', 'pca' or 'peer'. Methods to remove sample covariation in OUTRIDER.                                               ``autoencoder``
zScoreCutoff                  numeric    A non-negative number. Z scores (in absolute value) greater than this cutoff are considered as outliers.                               ``0``
padjCutoff                    numeric    A number between (0, 1] indicating the maximum FDR an event can have in order to be considered an outlier.                             ``0.05``
maxTestedDimensionProportion  numeric    An integer that controls the maximum value that the encoding dimension can take. Refer to the advanced options below.                  ``3``
============================  =========  =====================================================================================================================================  ======

Aberrant splicing dictionary
++++++++++++++++++++++++++++

============================  =========  ============================================================================================  ======
Parameter                     Type       Description                                                                                   Default/Examples
============================  =========  ============================================================================================  ======
groups                        list       Same as in aberrant expression.                                                               ``# see aberrant expression example``
minIds                        numeric    Same as in aberrant expression.                                                               ``1``
recount                       boolean    If true, it forces samples to be recounted.                                                   ``false``
longRead                      boolean    Set to true only if counting Nanopore or PacBio long reads.                                   ``false``
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
qcVcf                  character  Full path to the vcf file used for VCF-BAM matching                                                                       ``/path/to/qc_vcf.vcf.gz``
qcGroups               list       Same as “groups”, but for the VCF-BAM matching                                                                            ``# see aberrant expression example``
=====================  =========  ========================================================================================================================  ======


Creating the Sample Annotation Table
------------------------------------

For details on how to generate the sample annotation, please refer to the DROP manuscript. 
Here we provide some examples on how to deal with certain situations. For simplicity, we
do not include the other compulsory columns ``PAIRED_END``, ``COUNT_MODE``,
``COUNT_OVERLAPS`` and ``STRAND``.

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


Advanced options
----------------

A local copy of DROP can be edited and modified for uncovering potential issues or increasing outputs.
For example, the user might want to add new plots to the ``Summary`` scripts, or add
additional columns to the results tables.
Specifically, the number of threads allowed for a computational step can be modified by the user.

.. note::

    DROP needs to be installed from a local directory :ref:`otherversions` using ``pip install -e <path-to-drop-repo>``
    so that any changes in the code will be available in the next pipeline run.
    Any changes made to the R code need to be updated with ``drop update`` in the project directory.

The aberrant expression and splicing modules use a denoising autoencoder to
correct for sample covariation. This process reduces the fitting space to a 
dimension smaller than the number of samples N. The encoding dimension is optimized.
We recommend the search space to be at most N/3 for the aberrant expression, 
and N/6 for the aberrant splicing case. Nevertheless, the user can specify the 
denominator with the parameter ``maxTestedDimensionProportion``.

In order to influence which fields of the BAM files are imported, the user can 
provide a ``scanBamParam`` object. This will affect how the files are counted in 
the aberrant expression and splicing modules. Refer to the function's 
`documentation <https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/ScanBamParam>`_ for details.





