#'---
#' title: Results of FRASER analysis
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}--{annotation}" / "07_results.Rds")`'
#'  params:
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'   - outputDir: '`sm cfg.getProcessedResultsDir() + "/aberrant_splicing/datasets/"`'
#'   - padjCutoff: '`sm cfg.AS.get("padjCutoff")`'
#'   - zScoreCutoff: '`sm cfg.AS.get("zScoreCutoff")`'
#'   - deltaPsiCutoff: '`sm cfg.AS.get("deltaPsiCutoff")`'
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'   - assemblyVersion: '`sm cfg.genome.getBSGenomeVersion()`'
#'  threads: 10
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - add_HPO_cols: '`sm str(projectDir / ".drop" / "helpers" / "add_HPO_cols.R")`'
#'   - fdsin: '`sm cfg.getProcessedDataDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                 "padjBetaBinomial_theta.h5"`'
#'   - txdb: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/txdb.db"`'
#'   - gene_name_mapping: '`sm cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/gene_name_mapping_{annotation}.tsv"`'
#'   - spliceTypeSetup: '`sm cfg.AS.getWorkdir() + "/spliceTypeConfig.R"`'
#'   - addAnnotation:  '`sm cfg.AS.getWorkdir() + "/fds_annotation.R"`'
#'   - blacklist_19: '`sm cfg.AS.getWorkdir() + "/resource/hg19-blacklist.v2.bed.gz"`'
#'  output:
#'   - resultTableJunc: '`sm cfg.getProcessedResultsDir() + 
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv"`'
#'   - resultTableGene: '`sm cfg.getProcessedResultsDir() + 
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv"`'
#'   - fds: '`sm cfg.getProcessedResultsDir() + 
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/fds-object.RDS"`'                          
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$setup, echo=FALSE)
source(snakemake@input$add_HPO_cols)
source(snakemake@input$spliceTypeSetup, echo=FALSE)
source(snakemake@input$addAnnotation)
#source(snakemake@input$blacklist_19)
library(AnnotationDbi)

opts_chunk$set(fig.width=12, fig.height=8)

annotation    <- snakemake@wildcards$annotation
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir
outputDir  <- snakemake@params$outputDir 
assemblyVersion <- snakemake@params$assemblyVersion
#blacklist_folder <- snakemake@params$blacklist

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds_input <- loadFraserDataSet(dir=workingDir, name=dataset)
fds <- saveFraserDataSet(fds_input, dir=outputDir, name = paste(dataset, annotation, sep = '--'), rewrite = TRUE)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)

# Read annotation and match the chr style
txdb <- loadDb(snakemake@input$txdb)
orgdb <- fread(snakemake@input$gene_name_mapping)

seqlevelsStyle(orgdb$seqnames) <- seqlevelsStyle(fds)
seqlevelsStyle(txdb) <- seqlevelsStyle(fds)

# Annotate the fds with gene names
fds <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgdb, feature = 'gene_name', 
                              featureName = 'hgnc_symbol', keytype = 'gene_id')

# Add the junction annotations to the fds
fds <- createFDSAnnotations(fds, txdb)
#print(rowRanges(fds,type="j"))

mcols(fds, type="ss")$annotatedJunction <-  "none"
#print(rowRanges(fds, type="ss"))

# Extract results per junction
res_junc <- results(fds,
                    padjCutoff=snakemake@params$padjCutoff,
                    zScoreCutoff=snakemake@params$zScoreCutoff,
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff,
                    additionalColumns="annotatedJunction")
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')
saveFraserDataSet(fds, dir=outputDir)
#print(res_junc_dt)


# Add features 
if(nrow(res_junc_dt) > 0){
  
  # number of samples per gene and variant  
  res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
  res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
  res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
  
  # add colData to the results
  res_junc_dt <- merge(res_junc_dt, as.data.table(colData(fds)), by = "sampleID")
  res_junc_dt[, c("bamFile", "pairedEnd") := NULL]
} else{
  warning("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset.")
}


# Aggregate results by gene
if(length(res_junc) > 0){
  res_genes_dt <- resultsByGenes(res_junc) %>% as.data.table
  res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
  res_genes_dt[, c("bamFile", "pairedEnd") := NULL]
  
  # add HPO overlap information
  sa <- fread(snakemake@config$sampleAnnotation)
  if(!is.null(sa$HPO_TERMS)){
    if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
      res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = snakemake@params$hpoFile)
    }
  }
} else res_genes_dt <- data.table()

# Add UTR labels
res_junc_dt <- addUTRLabels(res_junc_dt, txdb)
#print(blacklist_19)
#blacklist_file <- blacklist_folder + "hg19-blacklist.v2.bed.gz"
#res_junc_dt <- addBlacklistLabels(res_junc, blacklist_file)

# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene)
