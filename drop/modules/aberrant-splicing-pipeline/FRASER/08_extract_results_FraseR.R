#'---
#' title: Results of FRASER analysis
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}--{annotation}" / "08_results.Rds")`'
#'  params:
#'   - workingDir: '`sm cfg.getProcessedResultsDir() + "/aberrant_splicing/datasets/"`'
#'   - padjCutoff: '`sm cfg.AS.get("padjCutoff")`'
#'   - zScoreCutoff: '`sm cfg.AS.get("zScoreCutoff")`'
#'   - deltaPsiCutoff: '`sm cfg.AS.get("deltaPsiCutoff")`'
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'  threads: 10
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - add_HPO_cols: '`sm str(projectDir / ".drop" / "helpers" / "add_HPO_cols.R")`'
#'   - fdsin: '`sm expand(cfg.getProcessedResultsDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/" +
#'                  "padjBetaBinomial_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'   - txdb: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/txdb.db"`'
#'  output:
#'   - resultTableJunc: '`sm cfg.getProcessedResultsDir() +
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv"`'
#'   - resultTableGene: '`sm cfg.getProcessedResultsDir() +
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$setup, echo=FALSE)
source(snakemake@input$add_HPO_cols)
library(AnnotationDbi)

annotation    <- snakemake@wildcards$annotation
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# Extract results per junction
res_junc <- results(fds, psiType=psiTypes,
                    padjCutoff=snakemake@params$padjCutoff,
                    zScoreCutoff=snakemake@params$zScoreCutoff,
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff)
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')


# Add features
if(nrow(res_junc_dt) > 0){

    # number of samples per gene and variant
    res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
    res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
    
    # add colData to the results
    res_junc_dt <- merge(res_junc_dt, as.data.table(colData(fds)), by = "sampleID")
    res_junc_dt[, c("bamFile", "pairedEnd", "STRAND") := NULL]
} else{
    warning("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset.")
}

# Extract results by gene
res_gene <- results(fds, psiType=psiTypes,
                    aggregate=TRUE, collapse=FALSE,
                    padjCutoff=snakemake@params$padjCutoff,
                    zScoreCutoff=snakemake@params$zScoreCutoff,
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff)
res_genes_dt   <- as.data.table(res_gene)
print('Results per gene extracted')

if(length(res_gene) > 0){
    res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
    res_genes_dt[, c("bamFile", "pairedEnd", "STRAND") := NULL]

    # add HPO overlap information
    sa <- fread(snakemake@config$sampleAnnotation, 
                colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
    if(!is.null(sa$HPO_TERMS)){
        if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
            res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = snakemake@params$hpoFile)
        }
    }
} else{
    res_genes_dt <- data.table()
    warning("The aberrant splicing pipeline gave 0 gene-level results for the ", dataset, " dataset.")
}

# Annotate results with spliceEventType and blacklist region overlap
# load reference annotation
library(AnnotationDbi)
txdb <- loadDb(snakemake@input$txdb)
    
# annotate the type of splice event and UTR overlap
res_junc_dt <- annotateSpliceEventType(result=res_junc_dt, txdb=txdb, fds=fds)
res_genes_dt <- annotateSpliceEventType(result=res_genes_dt, txdb=txdb, fds=fds)
    
# set genome assembly version to load correct blacklist region BED file (hg19 or hg38)
assemblyVersion <- snakemake@config$assemblyVersion
    
# annotate overlap with blacklist regions
res_junc_dt <- flagBlacklistRegions(result=res_junc_dt, 
                                    assemblyVersion=assemblyVersion)
res_genes_dt <- flagBlacklistRegions(result=res_genes_dt, 
                                     assemblyVersion=assemblyVersion)


# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene)
