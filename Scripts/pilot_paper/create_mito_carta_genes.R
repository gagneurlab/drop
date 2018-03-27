#'---
#' title: MitoCarta genes
#' author: Vicente Yépez
#' wb:
#'   input: "/s/project/mitoMultiOmics/raw_data/gene_info/Human.MitoCarta2.0.xls"
#'   output: 
#' output: 
#'   html_document
#'---
#'

#+ echo=F
library(data.table)
library(DT)
library(readxl)

# Input
dir_gene_info <- '/s/project/mitoMultiOmics/raw_data/gene_info/'
file_excel_mito_carta <- file.path(dir_gene_info, 'Human.MitoCarta2.0.xls')
mito_carta <- read_excel(file_excel_mito_carta, sheet = "A Human MitoCarta2.0") %>% as.data.table
# Output
file_out_mito_carta <- file.path(dir_gene_info, 'mito_carta_genes.tsv')

# Important cols
cols <- c(
    'Symbol', 
    'EnsemblGeneID',
    'Description',
    'Synonyms',
    'hg19_Chromosome',
    'Tissues',
    'MCARTA2.0_score',
    'MCARTA2_FDR')

mito_carta <- mito_carta[, cols, with=F]

setnames(mito_carta, "Symbol", "HGNC_GENE_NAME")
setnames(mito_carta, "EnsemblGeneID", "ENSEMBL_ID")
setnames(mito_carta, "hg19_Chromosome", "Chromosome")
setnames(mito_carta, "MCARTA2.0_score", "MCARTA_SCORE")
setnames(mito_carta, "MCARTA2_FDR", "FDR")

names(mito_carta) <- toupper(names(mito_carta))

write.table(mito_carta, file = file_out_mito_carta, quote = F, sep = "\t", row.names = F)


