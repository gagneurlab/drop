#'---
#' title: Results of FRASER analysis
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "07_results.Rds")`'
#'  params:
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'   - padjCutoff: '`sm cfg.AS.get("padjCutoff")`'
#'   - zScoreCutoff: '`sm cfg.AS.get("zScoreCutoff")`'
#'   - deltaPsiCutoff: '`sm cfg.AS.get("deltaPsiCutoff")`'
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'  threads: 10
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - add_HPO_cols: '`sm str(projectDir / ".drop" / "helpers" / "add_HPO_cols.R")`'
#'   - fdsin: '`sm cfg.getProcessedDataDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                 "padjBetaBinomial_theta.h5"`'
#'  output:
#'   - resultTableJunc: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/results/{dataset}_results_per_junction.tsv"`'
#'   - resultTableGene: '`sm cfg.getProcessedDataDir() + 
#'                          "/aberrant_splicing/results/{dataset}_results.tsv"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$setup, echo=FALSE)
source(snakemake@input$add_HPO_cols)

opts_chunk$set(fig.width=12, fig.height=8)

dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load data and annotate ranges with gene names
fds <- loadFraserDataSet(dir=workingDir, name=dataset)
GRCh <- ifelse(snakemake@config$genomeAssembly == 'hg19', 37, 
               ifelse(snakemake@config$genomeAssembly == 'hg38', 38,
                      error('Genome assembly must be either hg19 or hg38')))
fds <- annotateRanges(fds, GRCh = GRCh)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)

# Extract results per junction
res_junc <- results(fds,
                 padjCutoff=snakemake@params$padjCutoff,
                 zScoreCutoff=snakemake@params$zScoreCutoff,
                 deltaPsiCutoff=snakemake@params$deltaPsiCutoff,
                 additionalColumns=c("other_hgnc_symbol"))
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')
saveFraserDataSet(fds)

correctRes <- function(RT){
  rt <- copy(RT)
  rt[, bamFile := NULL]
  rt[, pairedEnd := NULL]
  rt[hgncSymbol == 'APOA1BP', hgncSymbol := 'NAXE']  # new gene name
  rt[hgncSymbol == 'C10orf2', hgncSymbol := 'TWNK']  # new gene name
  rt[hgncSymbol == 'CARKD', hgncSymbol := 'NAXD']  # new gene name
  return(rt)
}

# Add features 
if(nrow(res_junc_dt) > 0){
  
  # number of samples per gene and variant  
  res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
  res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
  res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
  
  # add colData to the results
  res_junc_dt <- merge(res_junc_dt, as.data.table(colData(fds)), by = "sampleID")
  res_junc_dt <- correctRes(res_junc_dt)
} else{
  warning("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset.")
}

# Aggregate results by gene
if(length(res_junc) > 0){
  res_genes_dt <- resultsByGenes(res_junc) %>% as.data.table
  res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
  res_genes_dt[, bamFile := NULL]
  res_genes_dt[, pairedEnd := NULL]
  
  # add HPO overlap information
  sa <- fread(snakemake@config$sampleAnnotation)
  if(!is.null(sa$HPO_TERMS)){
    if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
      res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = snakemake@params$hpoFile)
    }
  }
} else res_genes_dt <- data.table()

# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene)

