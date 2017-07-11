# Call

## SCARED
- consistency of results scared v1 to scared v2 incl Felix results
- patient report --> meta package summarizing 3 RNA defect strategies
- scared results: 
    - tidy table --> SummarizedExperiment
    - getResults returns table
- MAE pkg later use scripts for now
- devtools load folder <=> library()
- function names
    - think twice about it
    - use established from SummarizedExperiment
    - use same name as FraseR


## sample annotation
- implement as R pkg, not for public
- automated checks
- ask Tom again
    - build mapping from scratch
    - then match IDs with files we have
- updates on table: 
    - change root/source file
    - person shared betw HHZ and TUM, e.g. 1 day/week in HHZ
    - regular rebuild, e.g. daily
    - combine with rebuild upon trigger elements


## patient report
- focus on Exome and RNA
- modular elements
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



