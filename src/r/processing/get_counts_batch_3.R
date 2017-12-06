suppressPackageStartupMessages(source("src/r/config.R"))
exons_en = readRDS("resources/exons_en.Rds")
bf <- substr(scan("resources/201711_nadel_rna_seq_samples.txt", what = character()), 3, 100)
samples <- vapply(strsplit(bf, "/"), "[", "", 1)

bamfiles <- BamFileList(file.path(RAWDIR, "helmholtz", bf), yieldSize=2000000)
names(bamfiles) <- samples

my_bpparam = register_bplapply_for_clustering(slurm = F,
                                              workers = 40,
                                              threads = 2,
                                              memory = 8000,
                                              jobname = "count_rna"
)

starttime= Sys.time()

library(GenomicAlignments)
se_strand = summarizeOverlaps(
    exons_en,
    bamfiles,
    mode = 'IntersectionStrict',
    singleEnd = FALSE,
    ignore.strand = F, # TRUE,
    inter.feature = TRUE,     # TRUE, reads mapping to multiple features are dropped
    fragments = FALSE, 
    BPPARAM = my_bpparam    # SerialParam()
)
message('Processed all fibros in: ', format(Sys.time()- starttime))

saveRDS(se_strand, "/s/project/genetic_diagnosis/processed_data/Rds/se_batch3.Rds")