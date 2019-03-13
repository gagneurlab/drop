
vcfFile <- "/s/project/mitoMultiOmics/raw_data/helmholtz/74436/exomicout/paired-endout/processedData/vep_anno_74436_small_tmp.vcf.gz"
vcf <- readVcf(vcfFile)

source("Scripts/_functions/annotation_with_vep.R")

vep_param <- get_vep_params(version=94, num_forks=10, vcfFile=paste0(vcfFile, ".anno.vcf.gz"))
debug(ensemblVEP)
resCall <- ensemblVEP(vcfFile, vep_param)
resCall

vcfFileAnno <- "/s/project/mitoMultiOmics/raw_data/helmholtz/74436/exomicout/paired-endout/processedData/vep_anno_74436_small_tmp.vcf.gz.anno.vcf.gz"
vcf <- readVcf(vcfFileAnno)

utils::download.file("https://i12g-gagneurweb.in.tum.de/public/bugreports/bioc_variantAnnotation/example_no_anno.vcf.gz", "example_no_anno.vcf.gz")
utils::download.file("https://i12g-gagneurweb.in.tum.de/public/bugreports/bioc_variantAnnotation/example_vep_anno.vcf.gz", "example_vep.vcf.gz")

# plain vcf file
vcf <- readVcf("example_no_anno.vcf.gz")
colData(vcf)
dim(vcf)
str(seqlevels(vcf))

# annotated with VEP 
# contains very long line but no errors in the format
vcf <- readVcf("example_vep.vcf.gz")
colData(vcf)
dim(vcf)
str(seqlevels(vcf)[1])
str(seqlevels(vcf)[2])
str(seqlevels(vcf)[3])


# gr <- GRanges(seqnames = "chr2", ranges=IRanges(88999476-100, 89157282+100))
gr

ov <- findOverlaps(gr, vcf)
ov
to(ov)
rowRanges(vcf[to(ov)])
info(vcf[to(ov)])

# get commandline call
all_flags <- flags(vep_param)
all_flags <- all_flags[!sapply(all_flags, isFALSE)]
all_flags <- all_flags[!names(all_flags) %in% c('input_file', 'output_file', 'stats_file')]
config_string <- paste(names(all_flags), all_flags, sep = ' ', collapse = '\n')
# config_string <- gsub('TRUE|FALSE', '', config_string)
write(config_string, '~/Downloads/vep.config')

