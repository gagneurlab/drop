#'---
#' title: Exome candidates by variant prioritization
#' author: Daniel Bader
#' wb:
#'   input: [
#'     exome_candy: "/s/project/genetic_diagnosis/processed_results/variants_wes_candidates.RDS"
#'   ]
#'   output: 
#' output: 
#'   html_document
#'---
#'

#+ echo=F
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/diagnosis_tools/exome_candidates.R")
source("src/r/config.R")
file_exome_candy <- file.path(PROC_RESULTS, 'variants_wes_candidates.RDS')
file_disease_gene_anno <- file.path(RAWDIR, "gene_info/meta_disease_genes.tsv")

GENE_ANNO <- fread(file_disease_gene_anno, na.strings=c('NA',''))
exome_candy <- readRDS(file_exome_candy)


#+ echo=F
exome_candy <- exome_candy[order(FIBROBLAST_ID, chr, pos)]

#' compute size of variant
exome_candy[,ref:= as.character(ref)]
exome_candy[,alt:= as.character(alt)]
exome_candy[, var_size:=abs(nchar(ref)-nchar(alt))]

#+ echo=F
# subset columns
columns_to_show <- c(
    'sample',
    'hgncid',
    'gt',
    'mtype',
    'var_size',
    'mstype',
    'sift1',
    'pph1',
    'exacmaf',
    'sample_freq',
    'FIBROBLAST_ID',
    'chr',
    'pos',
    'ref',
    'alt'
)

# merge exome and disease gene info
exome_display_dt <- merge(
    GENE_ANNO[,.(HGNC_GENE_NAME, MIM_NUMBERS, DISEASE)], 
    exome_candy[, ..columns_to_show], 
    by.y='hgncid', 
    by.x='HGNC_GENE_NAME', 
    all.y=T
)

# round columns with numbers
columns_signif <- c('var_size', "sift1", "pph1", 'exacmaf', 'sample_freq')
for(j in columns_signif){
    exome_display_dt[, c(j):= list(signif(get(j), digits = 2))]
}

# renaming
setnames(exome_display_dt, 'mtype', 'variant type')
setnames(exome_display_dt, 'mstype', 'variant effect')
setnames(exome_display_dt, 'pph1', 'polyphen score')
setnames(exome_display_dt, 'sift1', 'sift score')
setnames(exome_display_dt, 'exacmaf', 'ExAC MAF')
setnames(exome_display_dt, 'sample', 'EXOME ID')




#' Variants are filtered for
#' 
#' - ExAC minor allele frequency < 0.001
#' - protein affecting
#' - potential bi-allelic effect on gene, i.e. homozygous or >=2 heterozygous variants in same gene
#' 

#+ echo=F
#head(exome_display_dt)
DT::datatable(
    exome_display_dt, 
    filter='top', 
    rownames = FALSE
    # ,options = list(scrollX = TRUE)
)



#+ END, echo=F


