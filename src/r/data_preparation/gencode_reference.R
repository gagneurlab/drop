# Download the latest GENCODE annotations file
# Author: vyepez
# It automatically recognizes the latest version, creates a folder (eg. gencode28)
# and exports the gff3 and fasta files
# On top, creates a gene annotation table with gene type and hgnc symbol


library(RCurl)
library(rtracklayer)
library(data.table)
library(magrittr)


# Variables ---------------------------------------------------------------

# 'human' or 'mouse'
organism <- "human/"

# output directory
out.dir <- "/s/genomes/human/hg19"


# Functions ---------------------------------------------------------------

ls_url <- function(url) {
    stopifnot(url.exists(url))
    out <- getURL(url, ftp.use.epsv = T, dirlistonly = TRUE)
    readLines(textConnection(out))
}


# Determine latest release ------------------------------------------------

ftp.url <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_", organism)
ftp.dirs <- ls_url(ftp.url)

ftp.dirs <- Filter(function(x) grepl("release", x), ftp.dirs)
releases <-sapply(strsplit(ftp.dirs, "_"), function(x) x[2])

ftp.dirs <- ftp.dirs[order(as.numeric(sub("\\D+", "", releases)), 
                           sub("\\d+", "", releases))]

last_release <- tail(ftp.dirs, 1)

message("Latest release is ", last_release)

# Make output directory
out.dir <- file.path(out.dir, paste0("gencode", strsplit(last_release, split = "_")[[1]][2]))
ifelse(!dir.exists(out.dir), dir.create(out.dir), FALSE)

# Download latest GFF3 and Fasta -----------------------------------------

ftp.url <- paste0(ftp.url, last_release, "/")
ftp.files <- ls_url(ftp.url)

gtf.file <- Filter(function(x) grepl("\\d\\.annotation.gff3", x), ftp.files)
gtf.url <- paste0(ftp.url, gtf.file)

fasta.file <- Filter(function(x) grepl("\\d\\.primary_assembly.genome.fa", x), ftp.files)
fasta.url <- paste0(ftp.url, fasta.file)

download.file(gtf.url, destfile = file.path(out.dir, gtf.file))
download.file(fasta.url, destfile = file.path(out.dir, fasta.file))


# Create gene annotation file -------------------------------------------
gtf <- readGFF(file.path(out.dir, gtf.file)) %>% as.data.table
gtf_gene <- gtf[type == "gene", .(seqid, source, type, start, end, strand, gene_id, gene_type, gene_name)]
write.table(gtf_gene, file.path(out.dir, "gene_annotation.tsv"), sep = "\t", quote = F, row.names = F)

