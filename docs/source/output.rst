Results and Output of DROP
===========================

DROP is intended to help researchers use RNA-Seq data in order to detect genes with aberrant expression,
aberrant splicing and mono-allelic expression. By simplifying the workflow process we hope to provide
easy to read and interpret html files and output files. This section is dedicated to explaining the relevant
results files. We will use the results of the ``demo`` to explain the files generated.::

    #install drop
    mamba create -n drop_env -c conda-forge -c bioconda drop
    conda activate drop_env
    
    mkdir drop_demo
    cd drop_demo
    drop demo
    
    snakemake -c1

Aberrant Expression
+++++++++++++++++++

html file
#########
Looking at the resulting ``Output/html/drop_demo_index.html`` we can see the ``AberrantExpression`` 
tab at the top of the screen. Following that the Overview tab contains links to the:  

* Counting Summaries 
    * For each aberrant expression group
        * split of local vs external sample counts
        * QC relating to reads and size factors for each sample
        * histograms relating to mean count distribution with different conditions
        * information about the expressed genes within each sample and as a dataset
* Outrider Summaries
    * For each aberrant expression group
        * the number of aberrantly expressed gene per sample
        * how batch correction is done and the resulting lack of batch effects
        * which samples contain outliers
        * results table
* Files
    * OUTRIDER files for each aberrant expression group
        * For each of these files you can follow the `OUTRIDER vignette for individual analysis <https://www.bioconductor.org/packages/devel/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf>`_. 
    * results.tsv files
        * For each aberrant expression group
            * a tsv file that contains the sampleID, hgnc gene symbol, pvalue and adjusted pvalue, and subsequent analysis points
                * this tsv file contains only the genes and samples that meet the cutoffs defined in the ``config.yaml``
                for ``padjCutoff`` and ``zScoreCutoff``

Aberrant Splicing
+++++++++++++++++

html file
##########
Looking at the resulting ``Output/html/drop_demo_index.html`` we can see the ``AberrantSplicing`` 
tab at the top of the screen. Following that the Overview tab contains links to the:  

* Counting Summaries 
    * For each aberrant splicing group
        * split of local (from internal BAM files) vs external sample counts
        * split of local vs merged with external sample splicing/intron counts
        * comparison of local and external log mean counts
        * histograms relating to junction expression before and after filtering and variability
* FRASER Summaries
    * For each aberrant splicing group
        * the number of samples, introns, and splice sites 
        * how batch correction is done and the resulting lack of batch effects
        * result table
* Files
    * FRASER files for each aberrant splicing group
        * For each of these files you can follow the `FRASER vignette for individual analysis <https://www.bioconductor.org/packages/devel/bioc/vignettes/FRASER/inst/doc/FRASER.pdf>`_. 
    * results.tsv files
        * For each aberrant splicing group
            * results.tsv 
                * this tsv file contains only significant junctions that meet the cutoffs defined in the ``config.yaml`` they are aggregated at the gene level. Any sample/gene pair is represented by only the most significant junction.
            * results_per_junction.tsv 
                * this tsv file contains only significant junctions that meet the cutoffs defined in the ``config.yaml`` they are aggregated at the junction level. 


Mono-allelic Expression
+++++++++++++++++++++++

html file
##########
Looking at the resulting ``Output/html/drop_demo_index.html`` we can see the ``MonoallelicExpression`` 
tab at the top of the screen. Following that the Overview tab contains links to the:  

* Results
    * For each mae group
        * the number of samples, unique genes, and aberrant events
        * a cascade plot that shows additional filters
            * MAE for REF: the monoallelic expression favors the reference allele 
            * MAE for ALT: the monoallelic expression favors the alternative allele 
            * rare: if ``add_AF`` is set to true in ``config.yaml`` must meet minimum AF set by ``max_AF``. Additionally it must meet the inner-cohort frequency ``maxVarFreqCohort`` cutoff
        * histogram of inner cohort frequency
        * summary of cascade plots and results table
* Files
    * Allelic counts
        * a directory containing the allelic counts of heterozygous variants
    * Results data tables of each sample (.Rds)
        * Rds objects containing the full results table regardless of MAE status
    * Significant MAE results tables
        * For each mae group
            * a link to the results tsv file. Only contains MAE results for the alternative allele
* Quality Control
    * QC Overview
        * For each mae group QC checks for DNA/RNA matching
* Analyze Individual Results
    * An example analaysis that can be run using the Rds objects linked in the files subsection
    * performed on the first mae sample