#'---
#' title: "RNA Variant Calling Summary: `r paste(snakemake@wildcards$dataset, snakemake@wildcards$annotation, sep = '--')`"
#' author: nickhsmith
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "{dataset}" / "{annotation}_RVC_summary.Rds")`'
#'  input:
#'   - singleVCF: '`sm createSingleVCF() `'
#'   - annotatedVCF: '`sm os.path.join(
#'                        cfg.processedResultsDir,
#'                        "rnaVariantCalling/out/batch_vcfs", "{dataset}",
#'                        "{dataset}_{annotation}.annotated.vcf.gz")`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] +
#'               "/rnaVariantCalling/{dataset}--{annotation}_summary.html"`'
#'  type: noindex
#'---  

#+ echo=FALSE
library(data.table)
library(ggplot2)
library(VariantAnnotation)
library(tMAE)
library(dplyr)
library(GenomicScores)

saveRDS(snakemake, snakemake@log$snakemake)
####
# Helper functions
####
.clean_vars <- function(x){
  x[which(x == "1/1" | x == "1|1")] <- "1/1"
  x[which(x == "1/0" | x == "1|0" | x == "0/1" | x == "0|1")] <- "0/1"
  x[which(x == "./." | x == ".|." | x == "0/0 "| x == "0|0")] <- "0/0"
  x[grepl("QD|FS|SnpCluster",x)] <- "Seq_filter"
  
  return(x)
}

vcf <- VariantAnnotation::readVcf(snakemake@input$annotatedVCF) # read the batch vcf
canonical_chr <- c(paste0("chr",c(1:22,"X","Y","M")),1:22,"X","Y","MT")
vcf <- vcf[seqnames(vcf) %in% canonical_chr]

vcf_dt <- as.data.table(geno(vcf)$GT)
vcf_dt[,FILTER := vcf@fixed$FILTER]

dt <- as.data.table(lapply(vcf_dt, .clean_vars))

if (snakemake@config$rnaVariantCalling$addAF){
  #tMAE code for adding gnomad frequency
  max_af_cutoff <- snakemake@config$rnaVariantCalling$maxAF
  pops <- c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax')
  maf_dt <- add_gnomAD_AF(granges(vcf),genome_assembly = snakemake@config$genomeAssembly) %>% as.data.table()

  # Compute the MAX_AF based on all provided population columns
  # return -1 if only NAs are present (to avoid a warning)
  maf_dt$MAX_AF <- apply(maf_dt[, ..pops], 1, 
                    FUN=function(x){ max(x, -1, na.rm=TRUE) })

  # Replace Inf/-1 with NA
  maf_dt[is.infinite(MAX_AF) | MAX_AF == -1, MAX_AF := NA]
  
  res <- cbind(dt,maf_dt[,"MAX_AF"])
  
  # calculate variant frequency within the cohort
  max_var_freq_cutoff <- snakemake@config$rnaVariantCalling$maxVarFreqCohort
  res$cohortFreq <- apply(res[,1:ncol(vcf)],1,function(x){
    1-sum(x == "0/0")/length(x) })

  # Label variants based on their MAX_AF
  res[FILTER == "PASS",FILTER:= "PASS_common"]
  res[MAX_AF <= max_af_cutoff & cohortFreq <= max_var_freq_cutoff & FILTER == "PASS_common",FILTER:= "PASS_rare"]
  
  #Reorder FILTER
  res$FILTER <- factor(res$FILTER,levels = c("Seq_filter","Mask","minALT","Mask;minALT","PASS_common","PASS_rare"))

  # drop AF and cohortFreq for plotting
  res[,MAX_AF := NULL]
  res[,cohortFreq := NULL]

} else{
  res <- dt
  #Reorder FILTER
  res$FILTER <- factor(res$FILTER,levels = c("Seq_filter","Mask","minALT","Mask;minALT","PASS"))
}


# Plot all filters
res_plot <- melt(res,id.vars = "FILTER",value.name = "GT")[GT != "0/0",.N,by = c("FILTER","variable","GT")] 
ggplot(res_plot, aes(x = FILTER, y = N,col = GT)) +
       geom_boxplot() +
       geom_text(data = res_plot[,median(N),by=c("FILTER","GT")],
           mapping = aes(x=FILTER,y= V1,label = V1, vjust = -0.5),position = position_dodge(0.9),show.legend = F,size = 3.5) +
       ylab("Variants per sample") + scale_x_discrete(guide = guide_axis(n.dodge = 2))

# Split res
res_plot[grepl("PASS",FILTER),FILTER := "PASS"]
res_plot[!grepl("PASS",FILTER),FILTER := "FILTERED"]

# Plot only Pass/Fail split
ggplot(res_plot, aes(x = FILTER, y = N,col = GT)) +
       geom_boxplot() +
       geom_text(data = res_plot[,median(N),by=c("FILTER","GT")],
           mapping = aes(x=FILTER,y= V1,label = V1, vjust = -0.5),position = position_dodge(0.9),show.legend = F,size = 3.5) +
       ylab("Variants per sample") + scale_x_discrete(guide = guide_axis(n.dodge = 2))
