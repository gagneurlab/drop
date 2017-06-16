# R function
#
# author: baderda


# source all R scripts in folder
load_rscripts_from_folder <- function(
    dirpath, pattern= "*.R", recursive=TRUE, ...
    ){
    filelist = list.files(
        dirpath, full.names = T, pattern= pattern, recursive=recursive, ...
        )
    tmpnull= sapply(filelist, source)
}
