# Create a link to the latest mito_disorder_prokisch_mayr table
# 1. Save Hans's Excel file into /s/project/mitoMultiOmics/raw_data/gene_info
# 2. Save the first page as tsv or csv YYYYmm_mito_disorder_genes_prokisch_mayr.tsv
# 3. Create a link: ln -s '/s/project/mitoMultiOmics/raw_data/gene_info/201902_mito_disorder_genes_prokisch_mayr.tsv' '/s/project/mitoMultiOmics/raw_data/gene_info/latest_mito_disorder_genes_prokish_mayr.tsv'
# 4. Run this script, hopefully there won't be any errors

require(data.table)

create_clean_prokisch_mayr_table <- function(input_file = NULL, output_file = NULL){

    dir_gene_info <- "/s/project/mitoMultiOmics/raw_data/gene_info/"
    if(is.null(input_file)){input_file <- file.path(dir_gene_info, 'latest_mito_disorder_genes_prokish_mayr.tsv')}
    if(is.null(output_file)){output_file <- file.path(dir_gene_info, 'mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv')}
    
    # process
    prokisch_mayr_dt <- as.data.table(read.delim(input_file, na.strings = '', stringsAsFactors = F))
    prokisch_mayr_dt <- clean_prokisch_mayr_gene_table(prokisch_mayr_dt)
    
    # rename
    setnames(prokisch_mayr_dt, gsub(" ", "_", toupper(names(prokisch_mayr_dt))))
    setnames(prokisch_mayr_dt, 'GENE', 'HGNC_GENE_NAME')
    setnames(prokisch_mayr_dt, 'ENERGY_METAB', 'DISEASE')
    setnames(prokisch_mayr_dt, 'ENTREZID', 'ENTREZ_GENE_ID')
    setnames(prokisch_mayr_dt, 'OMIM', 'MIM_NUMBER')

    # make list column to character
    prokisch_mayr_dt[, MIM_NUMBERS := paste(unlist(MIM_NUMBER),collapse = ','), by = 1:nrow(prokisch_mayr_dt)]
    prokisch_mayr_dt[, MIM_NUMBER := NULL]
    
    # fix entries
    prokisch_mayr_dt[, HGNC_GENE_NAME := toupper(HGNC_GENE_NAME)]
    prokisch_mayr_dt <- prokisch_mayr_dt[, lapply(.SD, function(j) gsub(' +$', '', j))]
    prokisch_mayr_dt[, OTHER_DESIGNATIONS := gsub('Other Designations: ', '', OTHER_DESIGNATIONS)]
    prokisch_mayr_dt[, ANNOTATION := gsub('Annotation: ','',ANNOTATION)]
    prokisch_mayr_dt[, OTHER_ALIASES := gsub('Other Aliases: ','',OTHER_ALIASES)]
    # prokisch_mayr_dt[, lapply(.SD, class)]
    
    prokisch_mayr_dt[HGNC_GENE_NAME == "C19ORF70,QIL1", HGNC_GENE_NAME := "C19ORF70"]
    
    prokisch_mayr_dt[CATEGORY == "Unclear function", CATEGORY := "Unclear Function"]
    prokisch_mayr_dt[CATEGORY == "homeostasis", CATEGORY := "Homeostasis"]
    
    prokisch_mayr_dt <- prokisch_mayr_dt[! HGNC_GENE_NAME %in% c("REMOVEDSACS", "REMOVEDHSPA9", "REMOVEDRTN4IP1", "?TOP3A")]
    
    # write clean output file
    write.table(prokisch_mayr_dt, file = output_file, quote=FALSE, sep='\t', row.names = F)
}


clean_prokisch_mayr_gene_table <- function(input_data){
    library(data.table)
    data <- copy(input_data)
    # find empty columns
    empty_coln <- t(data[,lapply(.SD, function(j) all(is.na(j)))])
    empty_coln <- as.data.table(empty_coln, keep.rownames=T)
    
    # remove empty columns
    data[, c(empty_coln[V1==T, rn]) := NULL]
    
    # clean colnames
    colnames <- gsub("(\\.| )+$", "", colnames(data))
    colnames <- gsub("\\.+", " ", colnames)
    colnames <- gsub("_", " ", colnames)
    colnames <- gsub("\\s*1.*", "", colnames)
    colnames <- gsub("Clinical synopsis", "CS", colnames)
    
    setnames(data, colnames)
    
    # rename some columns
    setnames(data, "GENE", "Gene")
    setnames(data, "First reported", "Year")
    setnames(data, "Gene ID", "EntrezID")
    setnames(data, "Gene Localisation", "Locus")
    setnames(data, "Name A K A", "Full Gene Name")
    setnames(data, "Mit Energy Metab", "Energy_Metab")
    setnames(data, "Associated disease phenotype s", "Associated_disease_phenotypes")
    
    
    # remove new lines within a cell
    data[, Associated_disease_phenotypes := gsub("\\s*\n\\s*", "; ", Associated_disease_phenotypes, perl = T)]
    # change number into real factor
    try(data[, Leigh := ifelse(Leigh == 1, "yes", NA)], silent = TRUE)
    try(data[, Neuro := ifelse(Neuro == '1', "yes",
        ifelse(Neuro == 'N', "no",
            ifelse(Neuro == '?', "unclear", 
                ifelse(Neuro == 'T', "tissue specific",
                    NA))))] , silent = TRUE)
    
    # Told by Hans
    data[Energy_Metab == "1?", Energy_Metab := '2']
    data[, Energy_Metab := ifelse(Energy_Metab == '1', "MITO",
        ifelse(Energy_Metab == '2', "maybe",
            ifelse(Energy_Metab == '3', "other",
                NA)))]
    try(data[, Cardiomyopathy := ifelse(Cardiomyopathy == 1, "major",
                                    ifelse(Cardiomyopathy == 2, "minor", NA))], silent = TRUE)
    
    omim_cols <- c("OMIM", "MIM", "CS", "CS2", "CS3", "CS4", "CS5", "CS6")
    # merge omim IDs
    omim_list <- apply(data[, ..omim_cols],1,function(x){
        ids <- unique(gsub(" +$", "", gsub("MIM: |[#*?%]| new OMIM number", "", na.omit(x))))
        ids[ids == 'No_MIM' | ids == ''] <- NA
        if(length(ids) > 0)
            return(as.character(na.omit(ids)))
        return(NA)
    })
    data[, OMIM := omim_list]
    
    # readable Entrez ID
    data[, EntrezID := gsub("ID: ", "", EntrezID)]
    data[, Chromosome := gsub("Chromosome: .?.?; Location: ", "", Chromosome)]
    data[, Locus := gsub(" ", "", Locus)]
    data[, Gene := gsub(" ", "", Gene)]
    
    # remove omim columns
    data[, c("Official Symbol", "Chromosome", "Gene Name", "MIM", 
        "CS", "CS2", "CS3", "CS4", "CS5", "CS6") := NULL]
    
    # remove authors
    data[, c("Count", "Sperl", "Prokisch", "Zeviani", "Thorburn", "Falk", "Taylor",
             "Authors I", "Authors II") := NULL]
    
    # remove other useless columns
    data[, c("Group", "ID", "Full Gene Name", "Statistik") := NULL]
    
    data <- data[!is.na(Gene)]
    
    return(data)
}


# Create the table
# dir_gene_info <- "/s/project/mitoMultiOmics/raw_data/gene_info/"
# create_clean_prokisch_mayr_table(file.path(dir_gene_info, 'latest_mito_disorder_genes_prokish_mayr.tsv'), file.path(dir_gene_info, 'mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv'))
