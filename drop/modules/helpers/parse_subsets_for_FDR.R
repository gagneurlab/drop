
# Function to extract subsets of genes for limited FDR correction from sample 
# annotation
parse_subsets_for_FDR <- function(yaml_file, sampleIDs){
    
    # if no file specific in config, return NULL
    if(is.null(yaml_file) || yaml_file == ""){
        return(NULL)
    }
    
    # check if file exists
    if(!file.exists(yaml_file)){
        stop("Cannot parse gene subsets: file ", yaml_file, " does not exist.")
    }
    
    # read in subsets from yaml
    subsets <- read_yaml(yaml_file)
    
    # expand subsets that should be tested for all samples
    subsets_expanded <- list()
    for(subset_name in names(subsets)){
        sub_ls <- subsets[[subset_name]]
        if(is.null(names(sub_ls))){
            # expand list to all samples if no sampleID specified in name
            sub_ls <- rep(list(unlist(sub_ls)), length(sampleIDs))
            names(sub_ls) <- sampleIDs
        } else{
            # remove sampleIDs not present in the current group
            sub_sids <- names(sub_ls)
            sub_ls <- sub_ls[sub_sids %in% sampleIDs]
        }
        # unless no sampleIDs present in this group, set to final subset list
        if(length(sub_ls) > 0){
            subsets_expanded[[subset_name]] <- sub_ls
        }
    }
    
    return(subsets_expanded)
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
