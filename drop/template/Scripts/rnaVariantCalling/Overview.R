#'---
#' title: RNA Variant Calling
#' author: Nick Smith
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "Overview.Rds")`'
#'  input:
#'    - vcfFilesBatch: '`sm expand(cfg.getProcessedDataDir() +
#'                      "/rnaVariantCalling/out/all_samples_haplocaller/" + 
#'                      "{batchID}_all_samples.genotyped.filtered_clean.vcf.gz",
#'                  batchID=cfg.RVC.groups)`'
#'    - countFiles: '`sm expand(cfg.getProcessedDataDir() +
#'                      "/rnaVariantCalling/out/sample_haplocaller/" + 
#'                      "{sample}/{sample}_variant_counts.txt",
#'                  sample=cfg.RVC.batchIDs,min_alt=getMinAlt())`'
#'    - vcfFilesMasked: '`sm expand(cfg.getProcessedDataDir() +
#'                      "/rnaVariantCalling/out/sample_haplocaller/" + 
#'                      "{sample}/{sample}.genotyped.filtered.basic{min_alt}.masked.vcf.gz",
#'                  sample=cfg.RVC.batchIDs,min_alt=getMinAlt())`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
library(data.table)
library(ggplot2)

saveRDS(snakemake, snakemake@log$snakemake)

#+ eval=TRUE, echo=FALSE
plot_table <- lapply(snakemake@input$countFiles,
  function(x){
  tmp_variants <- fread(x,sep = ",")
  }
)

plot_table <- rbindlist(plot_table)
plot_table$GT <- as.factor(plot_table$GT)
ggplot(melt(plot_table,measure.vars = 3:5),aes(x = variable, y= value,fill=GT)) + geom_boxplot() + 
  ylab("Number of Variants") + xlab("")

DT::datatable(plot_table, filter = 'top')
