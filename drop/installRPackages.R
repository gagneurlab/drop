options(timeout = 3600) # 1 hour for downloads

options(repos=structure(c(CRAN="https://cloud.r-project.org")), warn = -1)

if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
    BiocManager::install("remotes")
}
if (!requireNamespace('data.table', quietly = TRUE)) {
    install.packages('data.table')
}

suppressPackageStartupMessages(library(data.table))

# do not turn wanrings into errors. E.g. "Package XXX build for R 4.0.X"
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

args <- commandArgs(trailingOnly=TRUE)
if (file.exists(args[1])){
    packages <- fread(args[1], fill = TRUE)
} else {
    packages <- data.table(
        package=gsub("=.*", "", unlist(args)),
        version=gsub(".*=", "", unlist(args)),
        ref="")
    packages[package == version, c("min_version", "max_version") := ""]
}
installed <- as.data.table(installed.packages())

for (pckg_name in packages$package) {
    package_dt <- packages[package == pckg_name]
    pckg_name <- gsub(".*/", "", pckg_name)
    min_version <- package_dt$min_version
    max_version <- package_dt$max_version
    branch <- package_dt$ref
    
    if (!pckg_name %in% installed$Package || (min_version != "" && (compareVersion(
      installed[Package == pckg_name, Version], min_version) < 0 || compareVersion(
        installed[Package == pckg_name, Version], max_version) > 0))) {
        
        package <- package_dt$package
        message(paste("install", package))
        if(branch != ""){
            BiocManager::install(package, ask=FALSE, update=FALSE, ref=branch)
        } else{
            BiocManager::install(package, ask=FALSE, update=FALSE)    
        }
        message(paste("installed", package))
    }
}

options(warn = 0)
