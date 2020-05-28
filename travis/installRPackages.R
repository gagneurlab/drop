START_TIME <- Sys.time()

options(repos=structure(c(CRAN="https://cloud.r-project.org")))

if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
    install.packages('R.utils')
    BiocManager::install("remotes")
}


args <- commandArgs(trailingOnly=TRUE)
packages <- read.csv(args[1], stringsAsFactors = FALSE,
                     header = TRUE, sep = " ", comment.char = "#")
installed <- rownames(installed.packages())

install_packages <- function(packages) {
    for (i in 1:nrow(packages)) {
        
        package <- packages[i,1]
        pckg_name = tail(unlist(strsplit(package, split = "/")), n = 1)
        version <- packages[i, 'version']
        
        must_install <- TRUE
        
        # Do not install if already installed and version is at least the min required
        if (pckg_name %in% installed) {
            if((version == '' | compareVersion(as.character(packageVersion(pckg_name)), version) >= 0)){
                message(paste(pckg_name, "already installed"))
                must_install <- FALSE
            } 
        } 
        
        if(isTRUE(must_install)){
            message(paste("install", package))
            BiocManager::install(package)
            message(paste("installed", package))
        }
    }
}

maxTime <- max(30, (60*30 - difftime(Sys.time(), START_TIME, units="sec")))
R.utils::withTimeout(timeout=maxTime, {
    try({
        install_packages(packages)
    })
})

