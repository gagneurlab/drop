### Functions to add:
# 1. MITOCARTA T/F column
# 2. Hans mito list of disease associated genes
# 3. Disease associated genes from SAMPLE ANNOTATION
# 4. OMIM columns


# Read the latest gene annotation and subset it
# out.dir <- "/s/genomes/human/hg19/gencode28"   # This can change with a newer version
# gene_annot <- fread(file.path(out.dir, "gene_annotation.tsv"))
# gene_annot[, c("source", "type", "ID") := NULL]

require(data.table)
require(dplyr)
require(magrittr)

##############
### Add MITOCARTA T/F column
##############

# Read MitoCarta genes
add_mitocarta_col <- function(DT, gene_name_col = "gene_name"){
    DT <- copy(DT)
    dir_gene_info <- '/s/project/mitoMultiOmics/raw_data/gene_info/'
    mito_carta <- fread(file.path(dir_gene_info, 'mito_carta_genes.tsv'))
    
    # sds = setdiff(mito_carta$HGNC_GENE_NAME, gene_annot$gene_name) %>% sort
    
    # Create an alias table
    alias_dt_mc = data.table(v1 = c("ATP6", "ATP8","ND1","ND2","ND3","ND4","ND4L","ND5","ND6"))
    alias_dt_mc[, v2 := paste0("MT-", v1)]
    alias_dt_mc = rbind(alias_dt_mc, 
                    data.table(v1 = c("ACN9", "APOA1BP", "C10ORF2", "C6ORF57", "CARKD", "COA3", "COA4",
                                      "COX1", "COX2", "COX3", "CYTB", "HRSP12", "MTERF", "MTERFD1", "MTERFD2",
                                      "NRD1", "PET100", "PET112", "SLIRP", "SLMO1", "SLMO2", "TOMM70A", "XRCC6BP1",
                                      "C17ORF89", "ADCK3"), 
                              v2 = c("SDHAF3", "NAXE", "TWNK", "SDHAF4", "NAXD", "CCDC56", "CHCHD8", 
                                     "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "RIDA", "MTERF1", "MTERF3", "MTERF4",
                                     "NRDC", "C19ORF79", "GATB", "C14orf156", "PRELID3A", "PRELID3B", "TOMM70", "ATP23",
                                     "NDUFAF8", "COQ8A"))
                 )
    
    # Add MITOCARTA T/F column
    setnames(DT, gene_name_col, "gene_name")
    DT[, MITOCARTA := gene_name %in% c(mito_carta$HGNC_GENE_NAME, alias_dt_mc$v1, alias_dt_mc$v2)]
    setnames(DT, "gene_name", gene_name_col)
    
    return(DT)
}

# gene_annot <- add_mitocarta_col(gene_annot, gene_name_col = "gene_name")

##############
### Add Hans type of disease column
##############

# Read Hans gene disease list
add_hans_class <- function(DT, gene_name_col = "gene_name", return_all_info = TRUE){
    DT <- copy(DT)
    dir_gene_info <- '/s/project/mitoMultiOmics/raw_data/gene_info/'
    prokisch_mayr_dt <- fread(file.path(dir_gene_info, 'mitochondrial_disorder_genes_prokisch_mayr_cleaned.tsv'))
    
    setnames(prokisch_mayr_dt, "CATEGORY", "MITOGENE_CATEGORY")
    
    # Some of the genes have aliases, v1: must be the way they appear in Hans table
    alias_dt_hans = data.table(v1 = c("NAXE", "C19ORF70", "ATP5F1A", "ATP5F1E", "COQ8A", "COQ8B", "ATP5F1D", "PET100", "PARK2", "SPG20",
                                  "MRM2", "NDUFAF8", "RTN4IP1", "UQCC3", "TWNK", "FDX1L", "COA7", "USMG5", "ATP5MD_USMG5",
                                  "GATB", "NAXD"),
                               v2 = c("APOA1BP", "QIL1", "ATP5A1", "ATP5E", "ADCK3", "ADCK4", "ATP5D", "C19ORF79", "PRKN", "SPART",
                                  "FTSJ2", "C17ORF89", "NIMP", "C11ORF83", "C10ORF2", "FDX2",  "SELRC1", "ATP5MD", "ATP5MD",
                                  "CARKD", "PET112"))
    
    al = prokisch_mayr_dt[HGNC_GENE_NAME %in% alias_dt_hans$v1]
    al = merge(al, alias_dt_hans, by.x = "HGNC_GENE_NAME", by.y = "v1")
    al[, HGNC_GENE_NAME := v2]
    al[, v2 := NULL]
    
    pt = rbind(prokisch_mayr_dt, al)
    
    # sds = setdiff(prokisch_mayr_dt$HGNC_GENE_NAME, c(gene_annot$gene_name, alias_dt_hans$v1, alias_dt_hans$v2)) %>% sort
    
    setnames(pt, old = "DISEASE", new = "HANS_CLASS")
    
    # write.table(pt[HANS_CLASS == "MITO", sort(unique(HGNC_GENE_NAME))], "../mito-ncRNA/Data/mito_disease_genes.tsv", sep = "\n", row.names = F, col.names = F, quote = F)
    
    setnames(DT, gene_name_col, "gene_name")
    if(return_all_info == TRUE){
        DT <- left_join(DT, pt[,.(HGNC_GENE_NAME, HANS_CLASS, MITOGENE_CATEGORY, ASSOCIATED_DISEASE_PHENOTYPES)], 
                by = c("gene_name" = "HGNC_GENE_NAME")) %>% as.data.table 
    } else  {
        DT <- left_join(DT, pt[,.(HGNC_GENE_NAME, HANS_CLASS)], 
                        by = c("gene_name" = "HGNC_GENE_NAME")) %>% as.data.table
        DT[, MITO_DISEASE_GENE := FALSE]
        DT[HANS_CLASS == 'MITO', MITO_DISEASE_GENE := TRUE]
        DT[, HANS_CLASS := NULL]
    }
    
    setnames(DT, "gene_name", gene_name_col)
    
    return(DT)
}
# gene_annot = add_hans_class(gene_annot, gene_name_col = "gene_name")



##############
### Add disease genes columns
##############
add_disease_gene_info <- function(DT, gene_name_col = "gene_name"){
    dir_gene_info <- '/s/project/mitoMultiOmics/raw_data/gene_info/'
    disgene_dt <- fread(file.path(dir_gene_info, 'disease_genes_from_sample_anno.tsv'))
    
    # sds = setdiff(disgene_dt$HGNC_GENE_NAME, gene_annot$gene_name) %>% sort
    
    alias_dt_dis = data.table(v1 = c("APOA1BP", "PET100", "C10ORF2", "MNF1"),
                              v2 = c("NAXE", "C19ORF79", "TWNK", "UQCC2" ))
    
    al = disgene_dt[HGNC_GENE_NAME %in% alias_dt_dis$v1]
    al = merge(al, alias_dt_dis, by.x = "HGNC_GENE_NAME", by.y = "v1")
    al[, HGNC_GENE_NAME := v2]
    al[, v2 := NULL]
    
    pt = rbind(disgene_dt, al)
    
    setnames(pt, "DISEASE", "GENE_ASSO_DISEASE")
    
    setnames(DT, gene_name_col, "gene_name")
    DT <- left_join(DT, pt[,.(HGNC_GENE_NAME, GENE_ASSO_DISEASE)], 
                by = c("gene_name" = "HGNC_GENE_NAME")) %>% as.data.table
    setnames(DT, "gene_name", gene_name_col)
    
    return(DT)
}


##############
### Rahman's disease genes columns
##############
add_rahman_disease_col <- function(DT, gene_name_col = "gene_name"){
    rahman_table <- fread("/s/project/mitoMultiOmics/exomes1000/raw_data/genes_info/MitoDiseaseGenes(Rahman).csv")
    rahman_table[, Gene := toupper(Gene)]
    
    # v1 is the way it is in Rahman's table
    alias_dt = data.table(v1 = c("NAXE", "ATP5A1", "ATP5E", "COQ8A", "COQ8B", "ATP5D", "PET100", 
                                 "MRM2", "NDUFAF8", "RTN4IP1", "UQCC3", "TWNK", "FDX1L", "COA7"),
                          v2 = c("APOA1BP", "ATP5F1A", "ATP5F1E", "ADCK3", "ADCK4", "ATP5F1D", "C19ORF79", 
                                 "FTSJ2", "C17ORF89", "NIMP", "C11ORF83", "C10ORF2", "FDX2",  "SELRC1"))
    
    rt = rahman_table[Gene %in% alias_dt$v1]
    rt = merge(rt, alias_dt, by.x = "Gene", by.y = "v1")
    rt[, Gene := v2]
    rt[, v2 := NULL]
    
    rt = rbind(rahman_table, rt)
    
    setnames(DT, gene_name_col, "gene_name")
    DT[, rahman_disease_gene := gene_name %in% rt$Gene]
    setnames(DT, "gene_name", gene_name_col)
    
    return(DT)
}


##############
### Add OMIM columns
##############
# Adds 3 columns with OMIM number per gene, PINH: mode of inheritance, and PMIM: mim number of the phenotype
# Mind the script on "../mitomultiomics/src/r/functions/variant_handling/omim_parser.R"

add_omim_cols <- function(DT, gene_name_col = "gene_name", return_all_info = TRUE){
    omim_dt = readRDS("/s/project/mitoMultiOmics/db_data/omim-gene-pheno-cache.RDS")
    # omim_dt <- readRDS("../mitomultiomics/resource/omim_dt.Rds")
    omim_dt <- omim_dt[SYMBOL != "", .(SYMBOL, GMIM, PINH, PMIM)]
    omim_dt <- omim_dt[, SYMBOL := toupper(SYMBOL)]
    omim_dt <- omim_dt[, .SD[1], by = SYMBOL]   # In case of many GMIM, take the first only
    setnames(omim_dt, "GMIM", "OMIM")
    setnames(DT, gene_name_col, "gene_name")
    DT <- left_join(DT, omim_dt, by = c("gene_name" = "SYMBOL")) %>% as.data.table
    setnames(DT, "gene_name", gene_name_col)
    
    # return_all_info: 
    ## TRUE: returns all the 3 columns
    ## FALSE: returns a T/F OMIM col
    if(isFALSE(return_all_info)){
        DT[, OMIM_gene := !is.na(PMIM) & PMIM != ""]
        DT[, c("OMIM", "PINH", "PMIM") := NULL]
    } 
    return(DT)
}

# gene_annot = add_omim_cols(gene_annot, gene_name_col = "gene_name")

##############
### Respiratory chain complexes columns
##############
add_rcc_info <- function(DT, gene_name_col = "gene_name"){
    rcc_dt <- fread("/s/project/mitoMultiOmics/ncRNA/raw_data/ribosome_rcc_genes.tsv")
    rcc_dt = rcc_dt[Protein_Complex != "Ribosome"]
    rcc_dt[, Type := gsub("subunits", "subunit", Type)]
    rcc_dt[, Type := gsub("factors", "factor", Type)]
    rcc_dt[, RCC := paste(Protein_Complex, Type)]
    setnames(DT, gene_name_col, "gene_name")
    DT <- left_join(DT, rcc_dt[,.(Gene, RCC)], by = c("gene_name" = "Gene")) %>% as.data.table
    setnames(DT, "gene_name", gene_name_col)
    return(DT)
}

add_nuclear_mito_DNA <- function(DT, gene_name_col = "gene_name"){
    ge <- readRDS("../genetic_diagnosis/resources/gencode.v19_with_gene_name.Rds")
    ge[, nDNA_mtDNA := "nDNA"]
    ge[seqnames == "MT", nDNA_mtDNA := "mtDNA"]
    setnames(DT, gene_name_col, "gene_name")
    DT <- left_join(DT, ge[,.(gene_name, nDNA_mtDNA)], by = "gene_name") %>% as.data.table
    setnames(DT, "gene_name", gene_name_col)
    return(DT)
}

add_ensembl_id <- function(DT, gene_name_col = "gene_name", gene_id_col = "gene_id"){
    require(tidyr)
    ge <- readRDS("../genetic_diagnosis/resources/gencode.v19_with_gene_name.Rds")
    
    ge[, gene_name := toupper(gene_name)]
    
    # v1 the way it is in this annotation
    alias_dt = data.table(v1 = c("APOA1BP", "C19ORF70", "ATP5A1", "ATP5E", "ADCK3","ADCK4", "ATP5D", "PET100", 
                                 "FTSJ2", "NDUFAF8", "RTN4IP1", "UQCC3", "C10ORF2", "FDX1L", "SELRC1", "SEPN1", "PRUNE"),
                               v2 = c("NAXE", "QIL1", "ATP5F1A", "ATP5F1E", "COQ8A", "COQ8B", "ATP5F1D","C19ORF79", 
                                      "MRM2", "C17ORF89", "NIMP", "C11ORF83", "TWNK", "FDX2","COA7", "SELENON", "PRUNE1"))
    
    al = ge[gene_name %in% alias_dt$v1]
    al = merge(al, alias_dt, by.x = "gene_name", by.y = "v1")
    al[, gene_name := v2]
    al[, v2 := NULL]
    
    ge = rbind(ge, al)
    
    ge = separate(ge, gene_id, into = c("gene_id", "useless"), sep = "\\.", remove = T)
    setnames(DT, gene_name_col, "gene_name")
    DT <- left_join(DT, ge[,.(gene_name, gene_id)], by = "gene_name") %>% as.data.table
    
    # Add other ids
    DT[gene_name == "MT-ND4",  gene_id := "ENSG00000198886"]
    DT[gene_name == "MT-ATP6", gene_id := "ENSG00000198899"]
    DT[gene_name == "MT-TI",   gene_id := "ENSG00000210100"]
    setnames(DT, old = "gene_name", new = gene_name_col)
    setnames(DT, old = "gene_id", new = gene_id_col)
    
    return(DT)
}

add_gene_type <- function(DT, gene_name_col = "gene_name"){
    gt <- fread("/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv")
    setnames(DT, gene_name_col, "gene_name")
    DT <- left_join(DT, gt[,.(gene_name_unique, gene_type)], by = c("gene_name" = "gene_name_unique")) %>% as.data.table
    setnames(DT, "gene_name", gene_name_col)
    return(DT)
}

##############
### Add all columns
##############
add_all_gene_info <- function(DT, gene_name_col = "gene_name", mitocarta = T, hans = T, dis_genes = T, omim = T, rcc = F,
                              nDNA_mtDNA = F, gene_id = F, gene_id_col = "gene_id", gene_type = T, NA_as_dot = F){
    dt <- copy(DT)
    if(mitocarta == T) dt <- add_mitocarta_col(dt, gene_name_col)
    if(hans == T) dt <- add_hans_class(dt, gene_name_col)
    if(dis_genes == T) dt <- add_disease_gene_info(dt, gene_name_col)
    if(rcc == T) dt <- add_rcc_info(dt, gene_name_col)
    if(omim == T) dt <- add_omim_cols(dt, gene_name_col)
    if(nDNA_mtDNA == T) dt <- add_nuclear_mito_DNA(dt, gene_name_col)
    if(gene_id == T) dt <- add_ensembl_id(dt, gene_name_col, gene_id_col)
    if(gene_type == T) dt <- add_gene_type(dt, gene_name_col)
    if(NA_as_dot == T) dt[is.na(dt)] = "."
    return(dt)
}

# gene_annot <- add_all_gene_info(gene_annot, gene_name_col = "gene_name")
