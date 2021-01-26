options(repos=structure(c(CRAN="https://cloud.r-project.org")), warn = -1)
suppressPackageStartupMessages(library(data.table))

if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
    BiocManager::install("remotes")
}


args <- commandArgs(trailingOnly=TRUE)
packages <- fread(args[1], fill = TRUE)
packages <- packages[!startsWith(package, "#")]
installed <- as.data.table(installed.packages())

for (pckg_name in packages$package) {
    package_dt <- packages[package == pckg_name]
    pckg_name <- tail(unlist(strsplit(pckg_name, split = "/")), n = 1)
    version <- package_dt$version
    
    if (pckg_name %in% installed$Package &
      (version == "" || installed[Package == pckg_name, Version] == version)
    ) {
        #message(paste(pckg_name, "already installed"))
    } else {
        if (package_dt$bioconductor == TRUE) {
            INSTALL <- BiocManager::install
        } else {
            INSTALL <- install.packages
        }
        package <- package_dt$package
        message(paste("install", package))
        INSTALL(package)
        message(paste("installed", package))
    }
}

options(warn = 0)