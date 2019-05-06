#'---
#' title: Filter GTEx junction counts
#' author: vyepez
#' wb:
#'  input: 
#'  output:
#'  threads: 30
#'  type: script
#'---

saveRDS(snakemake, "tmp/filt_gtex_splicing.snakemake")
# snakemake <- readRDS("tmp/filt_gtex_splicing.snakemake")

suppressPackageStartupMessages({
  library(BiocParallel)
  library(SummarizedExperiment)
  devtools::load_all("../FraseR")
})

DIR_gtex_fraser <- "/s/project/gtex-processed/splicing_map"

gtex_dirs <- list.dirs(file.path(DIR_gtex_fraser, "savedObjects"), full.names = F)
if(any(grep("raw-", gtex_dirs))) gtex_dirs <- gtex_dirs[- grep("raw-", gtex_dirs)]
if(any(grep("-filtered", gtex_dirs))) gtex_dirs <- gtex_dirs[- grep("-filtered", gtex_dirs)]
gtex_dirs <- gtex_dirs[gtex_dirs != ""]

register(MulticoreParam(snakemake@threads))

bplapply(gtex_dirs, function(gd){
  message("loading ", gd)
  fds <- loadFraseRDataSet(DIR_gtex_fraser, gd)
  message("filtering ", gd)
  fds <- filterExpression(fds, filter = TRUE)
  message("saving ", gd)
  name(fds) <- paste(gd, 'filtered', sep = "-")
  # Important to save the whole dataset, not just the Rds file
  fds_new <- saveFraseRDataSet(fds, rewrite = TRUE)
})

