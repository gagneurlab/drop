# Call

## SCARED
- consistency of results scared v1 to scared v2 incl Felix results
    - get a SNP benchmark set to look for detection rate (Max)
- patient report --> meta package summarizing 3 RNA defect strategies
- scared results: 
    - tidy table --> SummarizedExperiment
    - getResults returns table
        - I would rather take results (R doesnt have the idea of setter/getter)
- MAE pkg later, use scripts for now
- devtools load folder <=> library()
    - document("./path/to/package/root")
    - install.packages("./path/to/package/root", repo=NULL, type="source")
    - load_all("./path/to/package/root")
- create your own class
    - let it inherit from SummariedExperiment
    - helps to implement easier specific functions for your package (addons)
- function names
    - think twice about it
    - use established from SummarizedExperiment/DESeq or other packages
    - use same name as FraseR
        - we could couple the dataset and inherit from another.
          then we do not need to implement some functions twice
            - samples 
            - samples <-
            - parallel
            - parallel <-
- setup Continuous Integration (gitlab-ci.yaml)
    - have a look at fraser/gitlab-ci.yaml
    - helps to see if package compiles, test succeeded, ...
    - tells you what is still missing for a final package :P 

## sample annotation
- implement as R pkg, not for public
- automated checks
    - contamination
    - sex
    - RNA/DNA matches
    - quality of reads, number of variants, ...
- ask Tom again
    - build mapping from scratch
        - based on IHG database
    - then match IDs with files we have
- updates on table: 
    - change root/source file
    - person shared betw HHZ and TUM, e.g. 1 day/week in HHZ
    - regular rebuild
        - eg: daily, snakemake (daily)
    - combine with rebuild upon trigger elements
        - issue, commits, ...
    - issue on gitlab for ID fix or other fixes in the data
- setup CI to get feedback on build status

## patient report
- focus on Exome and RNA
- modular elements
    - do not implement the visualization within the report
    - use the a object (widget) build by the packge to present the page
    - in server function call the provided server function from the widget object
- add OMIM genes to VIP genes: 
    - REST API
    - update every 6 months
    - mmo:src/r/functions/variant_handling/omim_parser.R


## sashimi plotly
- many singletons per junction
- reduce number of junctions plotted
- junctionSeq has quantitative splicing
- get Fraser example object to play: 
    - only exon junction counts
    - compute base pair counts dynamically


## Exome/ Genome Variants
- merge gvcf files with GATK --> one file for all
- work with VCF R object
- easy update with new samples



