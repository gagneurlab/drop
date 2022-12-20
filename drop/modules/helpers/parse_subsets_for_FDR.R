
# Function to extract subsets of genes for limited FDR correction from sample 
# annotation
parse_subsets_for_FDR <- function(sample_anno_file, sampleIDs,
                                    module=c("AE", "AS")){
    
    # read in sample annotation
    sa <- fread(sample_anno_file)
    sa <- sa[RNA_ID %in% sampleIDs]
    
    # define which column to extract from
    module <- match.arg(module)
    gene_subset_col <- paste0(module, "_GENES_TO_TEST")
    
    # check for presence of subset column
    if(!gene_subset_col %in% colnames(sa) || 
            all(sa[, get(gene_subset_col)] == "" |
                is.na(sa[, get(gene_subset_col)]))){
        return(NULL)
    }
    
    # extract sampleID and  file paths from sample annotation
    subset_files <- sa[, .(RNA_ID, get(gene_subset_col))]
    # split comma separated file into a separate rows (ignoring spaces if any)
    subset_files <- subset_files[, .(RNA_ID, unlist(tstrsplit(V2, ",(\\s)*")))]
    subset_files <- subset_files[!is.na(V2)]
    
    # read in gene subsets from files into list with respective setName
    subsets <- list()
    for(i in seq_len(subset_files[,.N])){
        sid <- subset_files[i, RNA_ID]
        file_path <- subset_files[i, V2]
        
        genesSubset <- fread(file_path)
        setName <- colnames(genesSubset)[1]
        genes <- genesSubset[,get(setName)]
        
        # throw error if no hashtag
        if(!grepl("^#", setName)){
            stop("Subset name missing in file: ", file_path)
        }
        # remove hashtag(s) and spaces from subset name
        setName <- gsub("^#(#)*(\\s)*", "", setName)
        
        if(setName %in% names(subsets)){
            prev_list <- subsets[[setName]]
            prev_list[[sid]] <- genes
            subsets[[setName]] <- prev_list
        } else{
            new_ls <- list(genes)
            names(new_ls) <- sid
            subsets[[setName]] <- new_ls
        }
    }
    
    # return subsets
    return(subsets)
}

# Convert element in subset gene lists from gene names to gene ids
convert_to_geneIDs <- function(subsets, gene_mapping_file){
    
    # if no subsets used, return NULL
    if(is.null(subsets)){
        return(NULL)
    }
    
    # otherwise map gene symbols to gene ids for OUTRIDER
    gene_annot_dt <- fread(gene_mapping_file)
    if(!is.null(gene_annot_dt$gene_name)){
        subsetsWithGeneIDs <- lapply(subsets, function(subset){
            geneIDs_mapped <- lapply(subset, function(genes){
                merge(data.table(gene_name=genes), 
                    gene_annot_dt[, .(gene_id, gene_name)],
                    by = 'gene_name', sort = FALSE)[, gene_id]
            })
        })
        return(subsetsWithGeneIDs)
    } else{
        warning("Cannot convert genes in subset to geneIDs, column gene_name ",
                "mising in annotation. Ignoring subsets for OUTRIDER.")
        return(NULL)
    }
    
}
