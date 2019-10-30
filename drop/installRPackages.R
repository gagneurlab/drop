args = commandArgs(trailingOnly=TRUE)
print(args)

if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
    BiocManager::install("remotes")
}

packages <- read.table(args[1], stringsAsFactors=FALSE)[,1]

for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        BiocManager::install(package)
        print(paste("installed", package))
    }
}
