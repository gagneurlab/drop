#'---
#' title: Collect all counts from FRASER Object
#' author: mumichae, vyepez, c-mertes
#' wb:
#'  log:
#'    - snakemake: '`sm str(tmp_dir / "AS" / "{dataset}" / "{genomeAssembly}--{annotation}_export.Rds")`'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '`sm cfg.getProcessedDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - fds: '`sm cfg.getProcessedDataDir() +
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'  output:
#'   - junction_counts: '`sm expand(cfg.exportCounts.getFilePattern(str_=True, expandStr=True) + "/raw_junction_counts.tsv.gz")`'
#'   - site_counts: '`sm expand(cfg.exportCounts.getFilePattern(str_=True, expandStr=True) + "/raw_site_counts.tsv.gz")`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

# 
# input
#
workingDir <- snakemake@params$workingDir
dataset    <- snakemake@wildcards$dataset

out_junction_counts <- snakemake@output$junction_counts
out_splice_site_counts <- snakemake@output$site_counts 


fds <- loadFraserDataSet(dir=workingDir, name=dataset)
chr_levels <- seqlevelsInUse(fds)

row_ranges <- as.data.table(rowRanges(fds))

junction_counts <- as.data.table(assays(fds)$rawCountsJ)
junction_counts <- cbind(junction_counts, row_ranges[, c("seqnames", "start", "end", "width", "strand", "startID", "endID")])
# Enforce order
junction_counts[, seqnames := factor(as.character(seqnames), levels = chr_levels)]
setorder(junction_counts, seqnames, start)

fwrite(junction_counts, file=out_junction_counts, quote=FALSE, row.names=FALSE, sep="\t", compress="gzip")

splice_sites <- rowData(nonSplicedReads(fds))
splice_site_counts <- as.data.table(assays(fds)$rawCountsSS)
splice_site_counts <- as.data.table(cbind(splice_sites, splice_site_counts))

splice_site_lookup <- unique(rbind(
  row_ranges[, .(spliceSiteID = startID, seqnames, position = start - 1)],  # donor positions
  row_ranges[, .(spliceSiteID = endID, seqnames, position = end)]      # acceptor positions
)) 

splice_site_lookup[, start := position]
splice_site_lookup[, end := start + 1]
splice_site_lookup[, width := 2]
splice_site_lookup[, position := NULL]

splice_site_counts <- merge(splice_site_counts, splice_site_lookup, by = "spliceSiteID", all.x = TRUE)        
splice_site_counts[, seqnames := factor(as.character(seqnames), levels = chr_levels)]
setorder(splice_site_counts, seqnames, start)



fwrite(splice_site_counts, file=out_splice_site_counts, quote=FALSE, row.names=FALSE, sep="\t", compress="gzip")

