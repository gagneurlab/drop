#'---
#' title: Expressed genes GTEx
#' author: vyepez
#' wb:
#'  input:
#'  output:
#'   - expressed_gtex: "/s/project/genetic_diagnosis/resource/gtex_expressed_genes.tsv"
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, "tmp/gtex_expressed.snakemake")
# snakemake <- readRDS("tmp/gtex_expressed.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(OUTRIDER)
  library(magrittr)
  library(dplyr)
})

# Read GTEx tissues
gtex_tissues <- fread("/s/project/genetic_diagnosis/resource/gtex_tissues.txt")

# Directory with ods objects of gtex tissues
gtex_dir <- "/s/project/scared/paper/revision/run0910-final/data/fitOutrider/NLas_TCNo"
gtex_files <- list.files(gtex_dir)

# Remove Kremer and simulations
gtex_files <- gtex_files[! gtex_files %in% c("Kremer_ODS.RDS", "SimulationNBinom_fitN_Q10_ODS.RDS", "SimulationNBinom_fitY_Q10_ODS.RDS", "SimulationNorm_fitN_Q10_ODS.RDS", "SimulationNorm_fitY_Q10_ODS.RDS")]

# Read ods and obtain expressed genes
DT <- lapply(gtex_files, function(gf){
  gtex_ods <- readRDS(file.path(gtex_dir, gf))
  dt <- data.table(Tissue_specific = gsub("_ODS.RDS", "", gf), gene = mcols(gtex_ods)$gene_symbol)
}) %>% rbindlist()

DT <- left_join(DT, gtex_tissues, by = 'Tissue_specific') %>% as.data.table()

fwrite(DT, snakemake@output$expressed_gtex)
