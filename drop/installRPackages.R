options(repos=structure(c(CRAN="https://cloud.r-project.org")))

if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
    BiocManager::install("remotes")
}

args <- commandArgs(trailingOnly=TRUE)
#package <- args[1]
packages <- read.table(args[1], stringsAsFactors=FALSE)[,1]

installed <- rownames(installed.packages())
for (package in packages) {
    # split package name from prefix
    pckg_name = tail(unlist(strsplit(package, split="/")), n=1)
    if (pckg_name %in% installed) {
        message(paste(pckg_name, "already installed"))
    } else {
        #BiocManager::install(package)
        message(paste("installed", package))
    }
}

