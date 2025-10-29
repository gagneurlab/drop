#'---
#' title: Collect all counts from FRASER Object
#' author: mumichae, vyepez, c-mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "{genomeAssembly}--{annotation}_export.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'

#'  input:
#'   - annotation: '`sm cfg.getProcessedDataDir() + "/preprocess/{annotation}/txdb.db"`'
#'   - splice_metrics: '`sm expand(cfg.getProcessedDataDir() +
#'                  "/aberrant_splicing/datasets/savedObjects/raw-local-{dataset}/{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'  output:
#'    - k_counts: '`sm expand(cfg.exportCounts.getFilePattern(str_=True, expandStr=True) + "/k_{metric}_counts.tsv.gz", metric=["j", "theta"])`'
#'    - n_counts: '`sm expand(cfg.exportCounts.getFilePattern(str_=True, expandStr=True) + "/n_{metric}_counts.tsv.gz", metric=["psi5", "psi3", "theta"])`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

library(AnnotationDbi)

# 
# input
#
annotation_file <- snakemake@input$annotation
workingDir <- paste(head(strsplit(snakemake@input$splice_metrics, "/")[[1]], -3), collapse = "/")
dataset    <- snakemake@wildcards$dataset

out_k_files <- snakemake@output$k_counts
out_n_files <- snakemake@output$n_counts

# Read annotation and extract known junctions
txdb <- loadDb(annotation_file)
introns <- unique(unlist(intronsByTranscript(txdb)))
introns <- keepStandardChromosomes(introns, pruning.mode = 'coarse')
length(introns)

# Read FRASER object, adapt chr style and subset to known junctions
fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset))
seqlevels(fds) <- seqlevelsInUse(fds)
seqlevelsStyle(fds) <- seqlevelsStyle(introns)[1]

# recalculate theta.h5 in Jaccard metric
if (fitMetrics(fds) == "jaccard"){
  fds <- calculatePSIValues(fds, type="theta")
} else {
  fds <- calculatePSIValues(fds, type="jaccard")
}

fds_known <- fds[unique(to(findOverlaps(introns, rowRanges(fds, type="j"), type="equal"))),]

# save k/n counts
sapply(c(out_k_files, out_n_files), function(i){
  ctsType <- toupper(strsplit(basename(i), "_")[[1]][1])
  psiType <- strsplit(basename(i), "_")[[1]][2]
  
  cts <- as.data.table(get(ctsType)(fds_known, type=psiType))
  grAnno <- rowRanges(fds_known, type=psiType)
  anno <- as.data.table(grAnno)
  anno <- anno[,.(seqnames, start, end, strand)]
  
  fwrite(cbind(anno, cts), file=i, quote=FALSE, row.names=FALSE, sep="\t", compress="gzip")
}) %>% invisible()
