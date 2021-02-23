#'---
#' title: Hyper parameter optimization
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "04_hyper.Rds")`'
#'  threads: 12
#'  input:
#'   - filter: '`sm expand(cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{{dataset}}/rawCounts{type}.h5",
#'                type=["J", "SS"]) `'
#'  output:
#'   - hyper: '`sm cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
suppressPackageStartupMessages({
  library(FRASER)
  library(tidyr)
})

if ("random_seed" %in% names(snakemake@config)){
  rseed <- snakemake@config$random_seed
  if(isTRUE(rseed)){
    set.seed(42)
  } else if (is.numeric(rseed)){
    set.seed(as.integer(rseed))
  }
}

#+ input
dataset  <- snakemake@wildcards$dataset
fds_file <- snakemake@input$filter[1]
BPPARAM  <- MulticoreParam(snakemake@threads)

# Set number of threads including for DelayedArray operations
register(BPPARAM)
DelayedArray::setAutoBPPARAM(BPPARAM)

# Force writing HDF5 files
options(FRASER.maxSamplesNoHDF5=-1)
options(FRASER.maxJunctionsNoHDF5=-1)

# Load PSI data
fds <- loadFraserDataSet(file=fds_file)

# Run hyper parameter optimization
implementation <- snakemake@config$aberrantSplicing$implementation
mp <- snakemake@config$aberrantSplicing$maxTestedDimensionProportion

# Get range for latent space dimension
a <- 2 
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

maxSteps <- 12
if(mp < 6){
  maxSteps <- 15
}

Nsteps <- min(maxSteps, b)
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique

for(type in psiTypes){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type, 
                            implementation=implementation,
                            q_param=pars_q,
                            plot = FALSE)
    fds <- saveFraserDataSet(fds)
}
fds <- saveFraserDataSet(fds)
file.create(snakemake@output$hyper)
