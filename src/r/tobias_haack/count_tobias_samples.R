suppressPackageStartupMessages(source("src/r/config.R"))
exons_op = readRDS("resources/exons_op.Rds")

source("src/r/tobias_haack/read_sample_anno.R")

samples <- sat[ASSAY == "RNASeq", ID_Links]
bf <- paste0(samples, "/RNAout/paired-endout/merged.rmdup.bam")

bamfiles <- BamFileList(file.path(RAWDIR, "helmholtz", bf), yieldSize=2000000)
names(bamfiles) <- samples

my_bpparam = register_bplapply_for_clustering(slurm = F, workers = 40, threads = 2, memory = 8000, jobname = "count_rna")


library(GenomicAlignments)

bplapply(1:length(bamfiles), function(i){
    se = summarizeOverlaps(
        exons_op,
        bamfiles[i],
        mode = 'IntersectionStrict',
        singleEnd = FALSE,
        ignore.strand = ! sat[ID_Links == names(bamfiles[i]), IS_RNA_SEQ_STRANDED],    # FALSE means strand specific
        inter.feature = TRUE,     # TRUE, reads mapping to multiple features are dropped
        fragments = FALSE   # SerialParam()
    )
    saveRDS(se, paste0("/s/project/genetic_diagnosis/processed_data/Rds/samples/", names(bamfiles[i]), ".Rds"))
},
BPPARAM = my_bpparam 
)
# message('Processed all fibros in: ', format(Sys.time()- starttime))
message('counting done')


# Test if the bam files exist
bplapply(1:length(bamfiles), function(i){
    file.exists(bamfiles[[i]]$path)
    
}, BPPARAM = my_bpparam
)
