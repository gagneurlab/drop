Preparing the Input Data
========================

Filling the Config File
-----------------------

The config file is in *YAML* format. It is composed of general and module-specific parameters. Variables are declared by writing the variable name followed by a colon, a space, and the value, for example:

``# full comment line``

``variable_1: /path/file.txt  # inline comment``

A variable can be of different types: string, numeric, list and dictionary. First, we describe each of those types of variables in *YAML* and then we provide a description and examples of each variable.


Types of variables in *YAML*
++++++++++++++++++++++++++++

==========  ===================================================================================================================================================================================================  ======================================
Type        Description                                                                                                                                                                                          Default/Examples
==========  ===================================================================================================================================================================================================  ======================================
boolean     A boolean is binary and can be true of false                                                                                                                                                         false
string      A text. Quotation marks are not needed, but we recommend them.                                                                                                                                       "/data/project1/sample_annotation.tsv"
numeric     An integer or a real number. A dot is used as a decimal separator.                                                                                                                                   0.05
list        A collection of elements of the same type. They are written with square brackets and its elements are separated by commas.                                                                           ['group1', 'group2']
dictionary  A dictionary (or named vector) containing key-value pairs. It is a collection of multiple elements where the key is a string and the value any type. Written in curly brackets and comma-separated.  {"f1": "file1.txt", "value": 0}
==========  ===================================================================================================================================================================================================  ======================================


Global parameters
+++++++++++++++++

================  ==========  =======================================================================================================================================  ==
Parameter         Type        Description                                                                                                                              Default/Examples
================  ==========  =======================================================================================================================================  ==
projectTitle      character   Title of the project to be displayed on the rendered HTML output                                                                         Project 1
htmlOutputPath    character   Full path of the folder where the HTML files are rendered                                                                                /data/project1/htmlOutput
root              character   Full path of the folder where the subdirectories processed_data and processed_results will be created containing DROP's output files     /data/project1
sampleAnnotation  character   Full path of the sample annotation table                                                                                                 /data/project1/sample_annotation.tsv
geneAnnotation    dictionary  A key-value list of the annotation name (key) and the full path to the GTF file (value). More than one annotation file can be provided.  {"anno1": /path/to/gtf.gtf}
tools             dictionary  A key-value list of different commands (key) and the command (value) to run them                                                         {gatkCmd: gatk, bcftoolsCmd: bcftools, samtoolsCmd: samtools}
================  ==========  =======================================================================================================================================  ==


Aberrant expression parameters
++++++++++++++++++++++++++++++

================  =======  =====================================================================================================================================  ==
Parameter         Type     Description                                                                                                                            Default/Examples
================  =======  =====================================================================================================================================  ==
minIds            numeric  Minimum number of samples required for the analysis.                                                                                   10
groups            list     DROP groups that should be executed in this module. If not specified, all groups are used.                                             ['group1', 'group2']
fpkmCutoff        numeric  A non-negative number indicating the minimum FPKM 5% of the samples per gene should have. If a gene has less it will be filtered out.  1 # suggested by OUTRIDER
zScoreCutoff      numeric  A non-negative number. Z scores (in absolute value) greater than this cutoff are considered as outliers.                               3
padjCutoff        numeric  A number between (0, 1] indicating the maximum FDR an event can have in order to be considered an outlier.                             0.05
useGeneNames      boolean  A boolean indicating whether to use hgnc gene symbols instead of ensembl gene ids.                                                     true
================  =======  =====================================================================================================================================  ==

Aberrant splicing parameters
++++++++++++++++++++++++++++

==============  =======  ============================================================================================  ==
Parameter       Type     Description                                                                                   Default/Examples
==============  =======  ============================================================================================  ==
minIds          numeric  Same as in aberrant expression.                                                               10
groups          list     Same as in aberrant expression.                                                               ['group1', 'group3']
deltaPsiCutoff  numeric  A non-negative number. Delta psi values greater than this cutoff are considered as outliers.  0.3 # suggested by FRASER
padjCutoff      numeric  Same as in aberrant expression.                                                               0.05
==============  =======  ============================================================================================  ==


Mono-allelic expression parameters
++++++++++++++++++++++++++++++++++

==================  ==========  ========================================================================================================================  ==
Parameter           Type        Description                                                                                                               Default/Examples
==================  ==========  ========================================================================================================================  ==
geneAssembly        character   Either hg19 or hg38, depending on the genome build used                                                                   hg19
fastaFile           character   Full path of a human reference fasta file                                                                                 /path/to/hg19.fa
groups              list        Same as in aberrant expression.                                                                                           ['group1', 'group2']
allelicRatioCutoff  numeric     A number between [0.5, 1) indicating the maximum allelic ratio allele1/(allele1+allele2) for the test to be significant.  0.8
padjCutoff          numeric     Same as in aberrant expression.                                                                                           0.05
maxAF               numeric     Maximum allele frequency (of the minor allele) cut-off. Variants with AF equal or below this number are considered rare.  0.001
addGnomAD           boolean     Whether or not to add the allele frequencies from gnomAD                                                                  true
qcVcf               dictionary  A key-value list of the chromosome format (key) and the full path to the VCF file used for matching samples (value)       {UCSC: /path/to/qc_ucsc.vcf.gz, NCBI: /path/to/qc_ncbi.vcf.gz}
qcGroups            list        Same as “groups”, but for the VCF-BAM matching                                                                            ['group1', 'group2']    
==================  ==========  ========================================================================================================================  ==


Creating the Sample Annotation Table
------------------------------------

For details on how to generate the sample annotation, please refer to the DROP paper. Here we will provide some examples.

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
