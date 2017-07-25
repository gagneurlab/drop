# R function
# 

create_clean_prokisch_mayr_table <- function(
    input_file=NULL, 
    output_file=NULL
){
    # set files
    # 
    dir_gene_info <- "/s/project/mitoMultiOmics/raw_data/gene_info/"
    if(is.null(input_file)){
        input_file <- file.path(
            dir_gene_info, 'mitochondrial_disorder_genes_prokisch_mayr.tsv'
        )
    }
    if(is.null(output_file)){
        output_file <- file.path(
            dir_gene_info, 
            'mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv'
        )
    }
    
    # process
    prokisch_mayr_dt <- as.data.table(read.delim(
        file = input_file, 
        na.strings = '', 
        stringsAsFactors = F
        ))
    prokisch_mayr_dt <- clean_prokisch_mayr_gene_table(prokisch_mayr_dt)
    
    # rename
    setnames(prokisch_mayr_dt, gsub(" ", "_", toupper(names(prokisch_mayr_dt))))
    setnames(prokisch_mayr_dt, 'GENE', 'HGNC_GENE_NAME')
    setnames(prokisch_mayr_dt, 'ENERGY_METAB', 'DISEASE')
    setnames(prokisch_mayr_dt, 'ENTREZID', 'ENTREZ_GENE_ID')
    setnames(prokisch_mayr_dt, 'OMIM', 'MIM_NUMBER')

    # make list column to character
    prokisch_mayr_dt[, MIM_NUMBERS:=paste(unlist(MIM_NUMBER),collapse = ','), by=ID]
    prokisch_mayr_dt[, MIM_NUMBER:=NULL]
    
    # fix entries
    # 
    prokisch_mayr_dt[, HGNC_GENE_NAME:=gsub(' ','',HGNC_GENE_NAME)]
    prokisch_mayr_dt <- prokisch_mayr_dt[,lapply(.SD, function(j) gsub(' +$', '', j))]
    # prokisch_mayr_dt[,lapply(.SD, class)]
    
    # write clean output file
    # 
    write_tsv(prokisch_mayr_dt, file = output_file)
}


clean_prokisch_mayr_gene_table <- function(input_data){
    library(data.table)
    data <- copy(input_data)
    # find empty columns
    empty_coln <- t(data[,lapply(.SD, function(j) all(is.na(j)))])
    empty_coln <- as.data.table(empty_coln, keep.rownames=T)
    empty_coln
    
    # remove empty columns
    data[,c(empty_coln[V1==T, rn]) := NULL]
    
    # clean colnames
    colnames <- gsub("(\\.| )+$", "", colnames(data))
    colnames <- gsub("\\.+", " ", colnames)
    colnames <- gsub("_", " ", colnames)
    colnames <- gsub("\\s*1.*", "", colnames)
    data <- setnames(data, colnames)
    
    # rename some columns
    setnames(data, "GENE", "Gene")
    setnames(data, "First reported", "Year")
    setnames(data, "Gene ID", "EntrezID")
    setnames(data, "Gene Localisation", "Locus")
    setnames(data, "Name A K A", "Full Gene Name")
    setnames(data, "Mit Energy Metab", "Energy Metab")
    
    # remove new lines within a cell
    data[,`Associated disease phenotype s`:=gsub("\\s*\n\\s*", "; ", `Associated disease phenotype s`, perl = T)]
    # change number into real factor
    data[,Leigh:=ifelse(Leigh == 1, "yes", NA)]
    data[,Neuro:=ifelse(Neuro == '1', "yes",
        ifelse(Neuro == 'N', "no",
            ifelse(Neuro == '?', "unclear", 
                ifelse(Neuro == 'T', "tissue specific",
                    NA
                ))))]
    data[,`Energy Metab`:=ifelse(`Energy Metab` == '1', "MITO",
        ifelse(`Energy Metab` == '2', "maybe",
            ifelse(`Energy Metab` == '3', "other",
                NA
            )))]
    
    # merge omim IDs
    omim_list <- apply(data[,.(OMIM,MIM,`Clinical synopsis`,CS2,CS3,CS4,CS5,CS6)],1,function(x){
        ids <- unique(gsub(" +$", "", gsub("MIM: |[#*?%]| new OMIM number", "", na.omit(x))))
        ids[ids == 'No_MIM' | ids == ''] <- NA
        if(length(ids) > 0)
            return(as.character(na.omit(ids)))
        return(NA)
    })
    data[,OMIM:=omim_list]
    
    # readable Entrez ID
    data[,EntrezID:=gsub("ID: ", "", EntrezID)]
    data[,Chromosome:=gsub("Chromosome: .?.?; Location: ", "", Chromosome)]
    data[,Locus:=gsub(" ", "", Locus)]
    
    # remove omim columns
    data[,c("Official Symbol", "Chromosome", "Gene Name", "MIM", 
        "Clinical synopsis", "CS2", "CS3", "CS4", "CS5", "CS6"):=NULL
        ]
    
    # sort it by ID
    setkey(data, ID)
    
    return(data)
}

