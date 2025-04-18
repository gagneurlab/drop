#'---
#' title: Estimating the optimal latent dimension
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "04_hyper.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 12
#'  input:
#'   - filter: '`sm expand(cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
#'  output:
#'   - hyper: '`sm expand(cfg.getProcessedDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

if ("random_seed" %in% names(snakemake@config)){
  rseed <- snakemake@config$random_seed
  if(isTRUE(rseed)){
    set.seed(42)
  } else if (is.numeric(rseed)){
    set.seed(as.integer(rseed))
  }
}

#+ input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load PSI data
fds <- loadFraserDataSet(dir=workingDir, name=dataset)
fitMetrics(fds) <- psiTypes

# Estimate optimal latent dimension using either OHT or grid search
implementation <- snakemake@config$aberrantSplicing$implementation
mp <- snakemake@config$aberrantSplicing$maxTestedDimensionProportion
oht <- snakemake@config$aberrantSplicing$useOHTtoObtainQ

if (isTRUE(oht)){
  message(date(), ": Using OHT implementation to determine optimal q ...")
  fds <- estimateBestQ(fds, type="jaccard",
                       useOHT=TRUE,
                       implementation=implementation,
                       plot = FALSE)
  metadata(fds)[["useOHTtoObtainQ"]] <- TRUE
} else{
  # Get range for latent space dimension
  a <- 2 
  b <- min(ncol(fds), nrow(fds)) / mp   # N/mp
  
  maxSteps <- 12
  if(mp < 6){
    maxSteps <- 15
  }
  
  Nsteps <- min(maxSteps, b)
  pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
  message(date(), ": Testing the following values of q to determine the optimal one: ",
          paste0(pars_q, collapse = ', '))
  
  for(type in psiTypes){
      message(date(), ": ", type)
      fds <- estimateBestQ(fds, type=type,
                           useOHT=FALSE,
                           implementation=implementation,
                           q_param=pars_q,
                           plot = FALSE)
      metadata(fds)[["useOHTtoObtainQ"]] <- FALSE
}}
fds <- saveFraserDataSet(fds)

# remove previous hyper.done files and create new one
outdir <- dirname(snakemake@output$hyper)
prevFilterFiles <- grep("hyper(.*)done", list.files(outdir), value=TRUE)
unlink(file.path(outdir, prevFilterFiles))
file.create(snakemake@output$hyper)
