#'---
#' title: Create datasets from annotation file
#' author: Christian Mertes, mumichae
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "00_defineDataset.Rds")`'
#'  params:
#'    - ids: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="RNA")`'
#'  input:
#'    - sampleAnnoFile: '`sm config["sampleAnnotation"]`'
#'  output:
#'    - colData: '`sm cfg.getProcessedDataDir() + 
#'                    "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'    - wBhtml:  '`sm config["htmlOutputPath"] + 
#'                    "/AberrantSplicing/annotations/{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages(library(data.table))

#+ input
outFile       <- snakemake@output$colData
annoFile      <- snakemake@input$sampleAnnoFile

#+ dataset name

name <- snakemake@wildcards$dataset
anno    <- fread(annoFile)

subset_ids <- snakemake@params$ids
annoSub <- anno[RNA_ID %in% subset_ids]
setnames(annoSub, "RNA_ID", "sampleID")
setnames(annoSub, "RNA_BAM_FILE", "bamFile")
setnames(annoSub, "PAIRED_END", "pairedEnd")
setcolorder(annoSub, unique(c("sampleID", "STRAND", "pairedEnd", "bamFile"), 
        colnames(annoSub)))

#'
#' ## Dataset: `r name`
#'
#+ echo=FALSE
finalTable <- annoSub

#'
#' ## Final sample table `r name`
#'
#+ savetable
DT::datatable(finalTable, options=list(scrollX=TRUE))

dim(finalTable)
write.table(x=finalTable, file=outFile, quote=FALSE, sep='\t', row.names=FALSE)
