#'---
#' title: Combine info on disease genes and solved patients
#' author: Daniel Bader
#' wb:
#'   input: '/s/project/mitoMultiOmics/raw_data/gene_info/disease_genes.tsv'
#'   output: "/s/project/mitoMultiOmics/raw_data/gene_info/disease_genes_from_sample_anno.tsv"
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

#'
#' ## File: list of genes by disease
#' 
disgene_dt <- fread(file.path(dir_gene_info, 'disease_genes.tsv'), na.strings = c("NA", ''))
#' * clean up disease names
setnames(disgene_dt, c('MITO', 'NBIA', 'other'))

#' * melt table
disgene_dt <- melt(disgene_dt, measure.vars=names(disgene_dt), variable.name='DISEASE', value.name='HGNC_GENE_NAME') %>% unique
disgene_dt <- disgene_dt[!is.na(HGNC_GENE_NAME)][order(HGNC_GENE_NAME)]

#' * Fix genes
disgene_dt <- disgene_dt[HGNC_GENE_NAME != "mtDNA"]
disgene_dt[, HGNC_GENE_NAME := toupper(HGNC_GENE_NAME)]
disgene_dt = disgene_dt[HGNC_GENE_NAME != "CLPB2"]  # CLPB2 doesn't exist in humans
disgene_dt = disgene_dt[HGNC_GENE_NAME != "TMEM126"]  # TMEM126 doesn't exist by itself. TMEM126A and TMEM126B exist.
disgene_dt = disgene_dt[HGNC_GENE_NAME != "TRMT1"]  # Confusion between TRMT1 and TRMU
disgene_dt = disgene_dt[HGNC_GENE_NAME != "IQSWC2"] # The gene doesn't exist

#' * New names
disgene_dt[HGNC_GENE_NAME == "C8ORF38", HGNC_GENE_NAME := "NDUFAF6"]
disgene_dt[HGNC_GENE_NAME == "TRX2", HGNC_GENE_NAME := "TXN2"]

#' * Fix diseases
disgene_dt[HGNC_GENE_NAME=='SPG7', DISEASE:='MITO']
disgene_dt[HGNC_GENE_NAME == "NBAS", DISEASE := "RALF"]
disgene_dt[HGNC_GENE_NAME %in% c("COASY", "SPG7"), DISEASE := "NBIA"]
disgene_dt[HGNC_GENE_NAME == "DNAJC3", DISEASE := "DNAJC3"]
disgene_dt[HGNC_GENE_NAME == "OCLN", DISEASE := "CNV"]
disgene_dt[HGNC_GENE_NAME %in% c("DHTKD1", "PRODH2"), DISEASE := "STOFFW"]
disgene_dt[HGNC_GENE_NAME == "SPG21", DISEASE := "NPC"]

#+
head(disgene_dt)
dim(disgene_dt)


library(stringr)
# Splits a column into multiple columns without having to tell how many
split_into_multiple <- function(column, pattern = ", ", into_prefix){
    cols <- str_split_fixed(column, pattern, n = Inf)
    # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
    cols[which(cols == "")] <- NA
    cols <- as.data.table(cols)
    names(cols) <- paste(into_prefix, 1:ncol(cols), sep = "_")
    return(cols)
}


sa <- SAMPLE_ANNOTATION[LAB == "PROKISCH"][!is.na(DISEASE) & !is.na(KNOWN_MUTATION)][, .(DISEASE, KNOWN_MUTATION)] %>% unique
sa_expanded = cbind(DISEASE = sa$DISEASE, split_into_multiple(sa$KNOWN_MUTATION, ";", "GENE"))
sa_melt = melt(sa_expanded, id.vars = "DISEASE", value.name = "HGNC_GENE_NAME")
sa_melt[, variable := NULL]
sa_melt <- na.omit(sa_melt)
sa_melt = sa_melt[DISEASE != "HEALTHY"]
barplot(table(sa_melt$DISEASE), las = 2)
dim(sa_melt)

dg <- rbind(disgene_dt, sa_melt) %>% unique

dg[, HGNC_GENE_NAME := toupper(HGNC_GENE_NAME)]
write_tsv(dg, file = file.path(dir_gene_info, "disease_genes_from_sample_anno.tsv"))


