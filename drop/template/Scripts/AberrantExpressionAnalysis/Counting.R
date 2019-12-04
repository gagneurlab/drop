#'---
#' title: RNA-Seq Counts
#' author: salazar, mumichae
#' wb:
#'  py:
#'   - |
#'    annotations = list(config["geneAnnotation"].keys())
#'    datasets = config["aberrantExpression"]["groups"]
#'  params:
#'    - tmpdir: '`sm drop.getTmpDir()`'
#'  input:
#'    - count_files: '`sm expand(parser.getProcDataDir() +
#'                    "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds",
#'                    annotation=annotations, dataset=datasets)`'
#'    - html: '`sm expand(config["htmlOutputPath"] +
#'              "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html",
#'              annotation=annotations, dataset=datasets)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "AberrantExpression_Counting.snakemake"))
# snakemake <- readRDS(".drop/tmp/AberrantExpression_Counting.snakemake")

groups <- names(snakemake@config$outrider_all)
anno_version <- names(snakemake@config$geneAnnotation)
html_file_dir <- gsub(snakemake@config$htmlOutputPath, '.', dirname(dirname(snakemake@input$html)))
summaries_titles <- sapply(anno_version, function(v) {
        paste0('[', groups ,'](', html_file_dir, '/', v, '/Summary_', groups, '.html){target="_blank"}', collapse = ' ')
    }
)

#' ## Count Summaries
#' `r paste0('Gene annotation version ', names(summaries_titles), ': ', summaries_titles, collapse = '\n')`
