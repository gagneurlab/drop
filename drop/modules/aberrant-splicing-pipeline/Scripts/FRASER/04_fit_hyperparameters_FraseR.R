#'---
#' title: Hyper parameter optimization
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 12
#'  input:
#'   - filter: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  output:
#'   - hyper: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_04.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_04.snakemake")

source("Scripts/_helpers/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load PSI data
fds <- loadFraserDataSet(dir=workingDir, name=dataset)

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

