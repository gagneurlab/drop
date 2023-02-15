#'---
#' title: Results of FRASER analysis
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}--{annotation}" / "08_results.Rds")`'
#'  params:
#'   - workingDir: '`sm cfg.getProcessedResultsDir() + "/aberrant_splicing/datasets/"`'
#'   - padjCutoff: '`sm cfg.AS.get("padjCutoff")`'
#'   - deltaPsiCutoff: '`sm cfg.AS.get("deltaPsiCutoff")`'
#'   - hpoFile: '`sm cfg.get("hpoFile")`'
#'   - ids: '`sm lambda w: sa.getIDsByGroup(w.dataset, assay="RNA")`'
#'   - genes_to_test: '`sm cfg.AS.get("genesToTest")`'
#'  threads: 10
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - add_HPO_cols: '`sm str(projectDir / ".drop" / "helpers" / "add_HPO_cols.R")`'
#'   - parse_subsets_for_FDR: '`sm str(projectDir / ".drop" / "helpers" / "parse_subsets_for_FDR.R")`'
#'   - sampleAnnoFile: '`sm config["sampleAnnotation"]`'
#'   - fdsin: '`sm expand(cfg.getProcessedResultsDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/" +
#'                  "padjBetaBinomial_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'   - txdb: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/txdb.db"`'
#'  output:
#'   - resultTableJunc: '`sm cfg.getProcessedResultsDir() +
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv"`'
#'   - resultTableGene_full: '`sm cfg.getProcessedResultsDir() +
#'                          "/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_gene_all.tsv"`'
#'   - resultTableGene_aberrant: '`sm cfg.getProcessedResultsDir() +
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

# read in subsets from sample anno if present (returns NULL if not present)
source(snakemake@input$parse_subsets_for_FDR)
fraser_sample_ids <- snakemake@params$ids
subsets <- parse_subsets_for_FDR(snakemake@params$genes_to_test,
                                 sampleIDs=fraser_sample_ids)

# Load fds and create a new one
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# Extract results per junction
res_junc <- results(fds, psiType=psiTypes,
                    padjCutoff=snakemake@params$padjCutoff,
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff,
                    subsets=subsets, 
                    fullSubset=FALSE)
res_junc_dt   <- as.data.table(res_junc)
colorder <- colnames(res_junc_dt[, !"FDR_set", with=FALSE])
res_junc_dt <- dcast(res_junc_dt, ... ~ FDR_set, value.var="padjust")
setnames(res_junc_dt, "transcriptome-wide", "padjust")
for(subset_name in names(subsets)){
    setnames(res_junc_dt, subset_name, paste0("padjust_", subset_name))
}
setcolorder(res_junc_dt, colorder)
print('Results per junction extracted')


# Add features
if(nrow(res_junc_dt) > 0){

    # number of samples per gene and variant
    res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
    res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
    
    # add colData to the results
    res_junc_dt <- merge(res_junc_dt, as.data.table(colData(fds)), by = "sampleID")
    res_junc_dt[, c("bamFile", "pairedEnd", "STRAND", "RNA_BAM_FILE", "DNA_VCF_FILE", "COUNT_MODE", "COUNT_OVERLAPS") := NULL]
} else{
    warning("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset.")
}

# Extract full results by gene
res_gene <- results(fds, psiType=psiTypes,
                    aggregate=TRUE, collapse=FALSE,
                    padjCutoff=NA, deltaPsiCutoff=NA, minCount=NA,
                    subsets=subsets, 
                    fullSubset=TRUE)
res_genes_dt   <- as.data.table(res_gene)
colorder <- colnames(res_genes_dt[, !c("padjust", "FDR_set"), with=FALSE])
res_genes_dt <- dcast(res_genes_dt[,!"padjust", with=FALSE], ... ~ FDR_set, value.var="padjustGene")
setnames(res_genes_dt, "transcriptome-wide", "padjustGene")
for(subset_name in names(subsets)){
    setnames(res_genes_dt, subset_name, paste0("padjustGene_", subset_name))
}
setcolorder(res_genes_dt, colorder)
print('Results per gene extracted')
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene_full)

# Subset gene results to aberrant
padj_cols <- grep("padjustGene", colnames(res_genes_dt), value=TRUE)
res_genes_dt <- res_genes_dt[do.call(pmin, c(res_genes_dt[,padj_cols, with=FALSE], 
                                                list(na.rm = TRUE))) <= snakemake@params$padjCutoff &
                                    abs(deltaPsi) >= snakemake@params$deltaPsiCutoff & 
                                    totalCounts >= 5,]

if(length(res_gene) > 0){
    res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
    res_genes_dt[, c("bamFile", "pairedEnd", "STRAND", "RNA_BAM_FILE", "DNA_VCF_FILE", "COUNT_MODE", "COUNT_OVERLAPS") := NULL]

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
res_junc_dt <- annotatePotentialImpact(result=res_junc_dt, txdb=txdb, fds=fds)
res_genes_dt <- annotatePotentialImpact(result=res_genes_dt, txdb=txdb, fds=fds)
    
# set genome assembly version to load correct blacklist region BED file (hg19 or hg38)
assemblyVersion <- snakemake@config$genomeAssembly
if(grepl("grch37", assemblyVersion, ignore.case=TRUE)){
    assemblyVersion <- "hg19"
}
if(grepl("grch38", assemblyVersion, ignore.case=TRUE)){
    assemblyVersion <- "hg38"
}
    
# annotate overlap with blacklist regions
if(assemblyVersion %in% c("hg19", "hg38")){
    res_junc_dt <- flagBlacklistRegions(result=res_junc_dt, 
                                        assemblyVersion=assemblyVersion)
    res_genes_dt <- flagBlacklistRegions(result=res_genes_dt, 
                                         assemblyVersion=assemblyVersion)
} else{
    message(date(), ": cannot annotate blacklist regions as no blacklist region\n", 
            "BED file is available for genome assembly version ", assemblyVersion, 
            " as part of FRASER.")
}

# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene_aberrant)
