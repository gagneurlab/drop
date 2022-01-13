#'---
#' title: RNA Variant Calling
#' author: Nick Smith
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "Overview.Rds")`'
#'  params:                                                                                                             
#'    - mafdb: '`sm cfg.genome.getMafDbName()`'
#'  input:
#'    - singleVCF: '`sm createSingleVCF() `'
#'    - annotatedVCF: '`sm expand(os.path.join(
#'                        cfg.processedResultsDir,
#'                        "rnaVariantCalling/out/batch_vcfs", "{dataset}",
#'                        "{dataset}_{annotation}.annotated.vcf.gz"), 
#'                    annotation = cfg.get("geneAnnotation"), dataset = cfg.RVC.groups) `'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
library(data.table)
library(ggplot2)
library(VariantAnnotation)
library(tMAE)
library(dplyr)

saveRDS(snakemake, snakemake@log$snakemake)
####
# Helper functions
####
.get_mafdb <- function(pkg_name){
  if(!requireNamespace(pkg_name, quietly=TRUE)){
    warning("The given MafDb is not installed: '", pkg_name, "'. We will do it now!")
    if(!requireNamespace("BiocManager", quietly=TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install(pkg_name, ask=FALSE)
  }
  
  mafdb <- getFromNamespace(pkg_name, pkg_name)
  mafdb
}
.clean_vars <- function(x){
  x[which(x == "1/1" | x == "1|1")] <- "1/1"
  x[which(x == "1/0" | x == "1|0" | x == "0/1" | x == "0|1")] <- "0/1"
  x[which(x == "./." | x == ".|." | x == "0/0 "| x == "0|0")] <- "0/0"
  x[grepl("QD|FS|SnpCluster",x)] <- "Seq_filter"
  
  return(x)
}


vcf <- VariantAnnotation::readVcf(snakemake@input$annotatedVCF[1]) # read the first batch
gr_vcf <- granges(vcf)

vcf_dt <- as.data.table(geno(vcf)$GT)
vcf_dt[,FILTER := vcf@fixed$FILTER]



dt <- as.data.table(lapply(vcf_dt, .clean_vars))

#tMAE code for adding gnomad frequency
mafdb <- .get_mafdb(snakemake@params$mafdb)
max_af_cutoff <- snakemake@config$mae$maxAF
populations <- c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax')

seqlevelsStyle(gr_vcf) <- seqlevelsStyle(mafdb)
chr_matching_index <- which(seqnames(gr_vcf) %in% seqnames(mafdb))
gr_vcf <- gr_vcf[chr_matching_index]

pt <- score(mafdb, gr_vcf, pop = populations) %>% as.data.table()
colnames(pt) <- populations
res <- cbind(dt[chr_matching_index], pt) %>% as.data.table()

# Compute the MAX_AF based on all provided population columns
# return -1 if only NAs are present (to avoid a warning)
res$MAX_AF <- apply(res[, ..populations], 1, 
                    FUN=function(x){ max(x, -1, na.rm=TRUE) })

# Replace Inf/-1 with NA
res[is.infinite(MAX_AF) | MAX_AF == -1, MAX_AF := NA]

# Label variants based on their MAX_AF
res[MAX_AF <= max_af_cutoff & FILTER == "PASS",FILTER:= "PASS_rare"]
res[MAX_AF > max_af_cutoff & FILTER == "PASS",FILTER:= "PASS_common"]
res[is.na(MAX_AF) & FILTER == "PASS",FILTER:= "PASS_rareNA"]

#Reorder FILTER
res$FILTER <- factor(res$FILTER,levels = c("Seq_filter","Mask","minALT","Mask;minALT","PASS_common","PASS_rareNA","PASS_rare"))

# drop populations and AF
res[,(populations) := NULL]
res[,MAX_AF := NULL]

# Plot all filters
ggplot(melt(res,id.vars = "FILTER")[value != "0/0",.N,by = c("FILTER","variable","value")], 
       aes(x = FILTER, y = N,col = value)) + geom_boxplot() +
      ylab("Variants per sample")


# Split res
res[grepl("PASS",FILTER),FILTER := "PASS"]
res[!grepl("PASS",FILTER),FILTER := "FILTERED"]

# Plot only Pass/Fail split
ggplot(melt(res,id.vars = "FILTER")[value != "0/0",.N,by = c("FILTER","variable","value")], 
       aes(x = FILTER, y = N,col = value)) + geom_boxplot() +
       ylab("Variants per sample")
