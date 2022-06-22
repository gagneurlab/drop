#'---
#' title: RNA Variant Calling
#' author: nickhsmith
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "RVC" / "Overview.Rds")`'
#'  params:
#'    - annotations: '`sm cfg.genome.getGeneVersions()`'
#'    - datasets: '`sm cfg.RVC.groups`'
#'    - htmlDir: '`sm config["htmlOutputPath"] + "/rnaVariantCalling"`'
#'  input:
#'    - functions: '`sm cfg.workDir / "Scripts/html_functions.R"`'
#'    - htmls:     '`sm expand(os.path.join(config["htmlOutputPath"],
#'                             "rnaVariantCalling",
#'                             "{annotation}/Summary_{dataset}.html"),
#'                          annotation = cfg.genome.getGeneVersions(), dataset = cfg.RVC.groups)`' 
#'    - annotated_vcfs: '`sm expand(os.path.join(
#'                            cfg.processedResultsDir,
#'                            "rnaVariantCalling/batch_vcfs", "{dataset}",
#'                            "{dataset}_{annotation}.annotated.vcf.gz"),
#'                          annotation = cfg.genome.getGeneVersions(), dataset = cfg.RVC.groups)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ include=FALSE
saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$functions)

#+ eval=TRUE, echo=FALSE
# get parameters
datasets <- sort(snakemake@params$datasets)
annotations <- snakemake@params$annotations
htmlDir <- snakemake@params$htmlDir

results_links <- sapply(
  annotations, function(x) build_link_list(
    file_paths = file.path(htmlDir, x, paste0('Summary_', datasets, '.html')),
    captions = datasets)
)


#'
#' **Datasets:** `r paste(datasets, collapse = ', ')`
#'
#' **Gene annotations:** `r paste(annotations, collapse = ', ')`
#'
#' ## RVC results
#' `r display_text(caption = 'Gene annotation_test:', links = results_links)`
#'
#' ## Files
#' * [single sample VCFs](`r file.path(snakemake@config$root, 'processed_results/rnaVariantCalling/sample_vcfs/')`)
#' * [batch annotated VCFs](`r file.path(snakemake@config$root, 'processed_results/rnaVariantCalling/batch_vcfs/')`)
