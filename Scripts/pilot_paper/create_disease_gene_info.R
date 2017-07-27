#'---
#' title: Combine info on disease genes and solved patients
#' author: Daniel Bader
#' wb:
#'   input: 
#'   output: 
#' output: 
#'   html_document:
#'     toc_float: yes
#'     toc: yes
#'---
#'



#+ echo=F
source("src/r/config.R")

#' 
#' # Collect disease genes
#' 
dir_gene_info <- '/s/project/mitoMultiOmics/raw_data/gene_info/'
dir_gene_info

file_meta_disease_gene <- file.path(dir_gene_info, 'meta_disease_genes.tsv')

#' 
#' ## File: list of genes by disease
#' 
disgene_dt <- fread(
    file.path(dir_gene_info, 'disease_genes.tsv'), 
    na.strings=c("NA", '')
)
#' * clean up disease names
setnames(disgene_dt, c('MITO', 'NBIA', 'other'))

#' * melt table
disgene_dt <- unique(melt(
    disgene_dt, 
    measure.vars=names(disgene_dt), 
    variable.name='DISEASE', 
    value.name='HGNC_GENE_NAME'
    ))
disgene_dt <- disgene_dt[!is.na(HGNC_GENE_NAME)][order(HGNC_GENE_NAME)]

#' * set SPG7 to MITO
disgene_dt[HGNC_GENE_NAME=='SPG7', DISEASE:='MITO']

#+
head(disgene_dt)
dim(disgene_dt)


#' 
#' ## File: paper disease gene list
#' 
paper_dt <- fread(
    file.path(dir_gene_info, 'kremer_bader_2016_biorxiv_mitochondrial_disease_genes.tsv'), 
    na.strings=c("NA", '')
)
setnames(paper_dt, c('HGNC_GENE_NAME', 'FULL_GENE_NAME', 'MIM_NUMBER', 'OMIM_LINK', 'ENTREZ_GENE_ID'))
paper_dt[is.na(MIM_NUMBER), OMIM_LINK:=NA]
paper_dt[, DISEASE:='MITO']
#+
head(paper_dt)
dim(paper_dt)


#' ### Merge
#' 
meta_disease_gene_dt <- merge(disgene_dt, paper_dt, all=T)

# remove spaces
meta_disease_gene_dt[, HGNC_GENE_NAME:=gsub(' ','',HGNC_GENE_NAME)]
# remove trailing spaces
meta_disease_gene_dt <- meta_disease_gene_dt[,lapply(.SD, function(j) gsub(' +$', '', j))]



#' 
#' ## File: manually curated list by Prokisch and Mayr
#' 
create_clean_prokisch_mayr_table()
#' 
prokisch_mayr_dt <- fread(
    file.path(dir_gene_info, 'mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv')
)
#+
sort(names(prokisch_mayr_dt))


#'
#' ### Merge
#'
dt <- merge(meta_disease_gene_dt, prokisch_mayr_dt, all=T, by='HGNC_GENE_NAME')
dim(dt)

#' * drop full gene names; can be inferred later, if needed
dt[,FULL_GENE_NAME.x:=NULL]
dt[,FULL_GENE_NAME.y:=NULL]
dt[,ID:=NULL]

#' * fill NAs in duplicated columns bi-directional
icols <- c('HGNC_GENE_NAME', 'DISEASE.x', 'DISEASE.y', 'ENTREZ_GENE_ID.x', 'ENTREZ_GENE_ID.y')
merged_cols <- gsub('\\.x', '', grep('\\.x', names(dt), value=T))

for(coln in merged_cols){
    dt[is.na(get(paste0(coln,'.x'))) | is.na(get(paste0(coln,'.y'))), 
        eval(coln) := na.omit(c(get(paste0(coln,'.x')), get(paste0(coln,'.y')))), 
        by=HGNC_GENE_NAME
        ]    
}


#' * check duplicated columns for mismatches
#'
#'   * DISEASE
#+ results='hold'
dt[DISEASE.x != DISEASE.y, icols, with=F]
dt[is.na(DISEASE), DISEASE:= DISEASE.x]
#' DISEASE --> keep coln.x

#'   * ENTREZ gene id
#+ results='hold'
dt[ENTREZ_GENE_ID.x != ENTREZ_GENE_ID.y, c(icols,merged_cols), with=F]
dt[is.na(ENTREZ_GENE_ID), ENTREZ_GENE_ID:= ENTREZ_GENE_ID.x]
#' ENTREZ_GENE_ID --> all fine
#' 
#' * remove duplicated columns
dupl_coln <- grep('\\.[xy]', names(dt), value=T)
dt[,c(dupl_coln):=NULL]


#' MIM numbers
#' 
#' * add single MIM_NUMBERS to list column, where list is NA
dt[is.na(MIM_NUMBERS) &  !is.na(MIM_NUMBER), MIM_NUMBERS:=MIM_NUMBER]
#' * check single number is always contained in list
#' 
# make list
dt[, list_mim:=strsplit(MIM_NUMBERS, ',')]
# annotate length
dt[, num_mim_numbers:=length(list_mim[[1]]), by=HGNC_GENE_NAME]
# compute matching
dt[, mimX_in_mimY:=MIM_NUMBER %in% unlist(list_mim), by=HGNC_GENE_NAME]
#
nrow(
    dt[!is.na(MIM_NUMBER) & mimX_in_mimY==F, 
        .(HGNC_GENE_NAME, MIM_NUMBER, list_mim)]
    )==0
#' --> all single MIM_NUMBERS are contained in the list column
#' 
#' * remove helper columns
dt[,c('MIM_NUMBER', 'list_mim', 'num_mim_numbers', 'mimX_in_mimY'):=NULL]


#' Final merged columns
#' 
#+ results='hold'
meta_disease_gene_dt <- unique(dt)
write_tsv(meta_disease_gene_dt, file=file_meta_disease_gene)
#
dim(meta_disease_gene_dt)
meta_disease_gene_dt[,.N, DISEASE]
sort(names(meta_disease_gene_dt))



#' 
#' ## FILE: mito associated genes
#' 

mito_asso_dt <- fread(
    file.path(dir_gene_info, 'mito_associated_genes.tsv'), 
    na.strings=c("NA", '')
)
mito_asso_dt <- mito_asso_dt[,lapply(.SD, function(j) gsub(' +$', '', j))]
mito_asso_dt[, c('alternativeName', 'encodingOrigin'):=NULL]
setnames(mito_asso_dt, 'GeneName', 'HGNC_GENE_NAME')
setnames(mito_asso_dt, 'Complex', 'FUNCTION')

#' Are there new gene to add?
#' 
nrow(mito_asso_dt[!HGNC_GENE_NAME %in% meta_disease_gene_dt$HGNC_GENE_NAME])>0

#' **==> GENE TABLE IS COMPLETE**

#'
#' # Consistency checks
#'
#' ## Validate gene names 
#' 
#' Are the duplicates?
genes_dupl <- meta_disease_gene_dt[,.N,HGNC_GENE_NAME][N>1, HGNC_GENE_NAME]
genes_dupl
meta_disease_gene_dt[HGNC_GENE_NAME %in% genes_dupl]
#' 
#' --> duplicate origins from Prokisch Mayr table
#' 

#' Load HGNC annotation
#' 
hgnc_dt <- fread(
    '/s/genomes/human/hg19/hgnc/gene_with_protein_product_20151027.txt', 
    na.strings=c("NA", '')
)

sort(setdiff(meta_disease_gene_dt$HGNC_GENE_NAME, hgnc_dt$`Approved Symbol`))


