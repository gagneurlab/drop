#'---
#' title: Create datasets from annotation file
#' author: Christian Mertes
#' wb:
#'  params:
#'    - ids: '`sm parser.fraser_ids`'
#'    - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'    - fileMappingFile: '`sm parser.getProcDataDir() + "/file_mapping.csv"`'
#'  input:
#'    - sampleAnnoFile: '`sm config["sampleAnnotation"]`'
#'  output:
#'    - colData: '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_00.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_00.snakemake")

#+ load main config, echo=FALSE
source("Scripts/_helpers/config.R", echo=FALSE)

#+ input
outFile       <- snakemake@output$colData
annoFile      <- snakemake@input$sampleAnnoFile
fileMapFile   <- snakemake@params$fileMapping

#+ dataset name

name <- snakemake@wildcards$dataset
anno    <- fread(annoFile)
mapping <- fread(fileMapFile)

subset_ids <- snakemake@params$ids[[name]]
annoSub <- anno[RNA_ID %in% subset_ids]
colData <- merge(
  annoSub[,.(sampleID = RNA_ID, STRAND, PAIRED_END)],
  mapping[FILE_TYPE == "RNA_BAM_FILE", .(sampleID=ID, bamFile=FILE_PATH)])
colData <- unique(colData)  # needed in case of DNA replicates

write_tsv(colData, file=outFile)
