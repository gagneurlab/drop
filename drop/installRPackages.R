options(repos=structure(c(CRAN="https://cloud.r-project.org")))

if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
    BiocManager::install("remotes")
}

args = commandArgs(trailingOnly=TRUE)
package = args[1]

if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
    message(paste("installed", package))
} else {
    message(paste(package, "already installed"))
}

