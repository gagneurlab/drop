# R functions
# author: mertes, baderda
#
############################################


#' List available nodes
#' 
get_slurm_nodes= function(as_vector=FALSE){
    
    # you need to be on a slurm node for that command
    sinfo_out= system("sinfo", intern=TRUE)
    
    # nodes are in coln 6, drop coln header with -1
    nodelist= sapply(strsplit(sinfo_out, ' +'), '[', 6)[-1]
    
    # output formatting
    if(as_vector){
        return(nodelist)
    }else{
        return(paste(nodelist, collapse=','))
    }
}

#'
#' register bplapply with slurm or multicore
#'
register_bplapply_for_clustering <- function(
    slurm = TRUE,
    workers = 16, 
    threads = 8,
    memory = 8000, 
    jobname = "mitoSlurmBP", 
    workdir = getwd(),
    time = "10:00:00"
# TODO		,nodelist= get_slurm_nodes()
){
    library(BatchJobs)
    library(BiocParallel)
    
    # check for file
    file_slurm_template <- "./slurm.tmpl"
    stopifnot(file.exists(file_slurm_template))
    
    # create output dir
    if(grepl("/", jobname)){
        jobname <- gsub(".*/", "", jobname)
    }
    log_out_dir <- file.path(DATADIR, paste0("tmp_", Sys.info()['user']), "slurm_output/", jobname)
    if(!file.exists(log_out_dir)){
        dir.create(log_out_dir, recursive = TRUE)
    }
    
    if(slurm){
        # register slurm 
        funs <- makeClusterFunctionsSLURM(file_slurm_template)
        
        # setConfig(debug = TRUE)
        setConfig(raise.warnings = TRUE)
        db.options = list(pragmas = c("busy_timeout=10000", "journal_mode=WAL"))
        setConfig(db.options = db.options)
        setConfig(fs.timeout = 100)
        
        # slurm param object
        param <- BatchJobsParam(
            workers, 
            cluster.functions=funs, 
            resources=list(
                ncpus=threads, 
                memory=sprintf("%d", memory), 
                wd = workdir,
                out=file.path(log_out_dir, "slurm_out_%J.out"),
                err=file.path(log_out_dir, "slurm_out_%J.err"),
                time = time
# TODO						, nodelist= nodelist
            ), 
            jobname = jobname,
            #work.dir = file.path(tempdir(), 'biocParallel')
            work.dir = workdir
        )
        
    } else {
        param <- MulticoreParam(workers * threads)
    }
    
    # register it
    register(param)
    return(param)
}

