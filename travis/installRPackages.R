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
    
    pckg_name = tail(unlist(strsplit(packages[i,1], split = "/")), n = 1)
    version <- packages[i, 'version']
    print(version)
    right_version <- (is.na(version) | compareVersion(as.character(packageVersion(pckg_name)), version) >= 0)
    print(right_version)
    # right_version <- TRUE
    
    if (pckg_name %in% installed & isTRUE(right_version)) {
        message(paste(pckg_name, "already installed"))
    } else {
        if (packages[i,2] == TRUE) {
            INSTALL <- BiocManager::install
        } else {
            INSTALL <- install.packages
        }
        package <- packages[i,1]
        message(paste("install", package))
        INSTALL(packages[i,1])
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

