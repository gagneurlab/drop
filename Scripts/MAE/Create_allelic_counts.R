devtools::load_all("../mae/")
pc <- fread('/s/project/genetic_diagnosis/processed_data/mae_files.tsv')

vcfs <- pc$exome_vcf_file
rnas <- pc$RNA_file

BPPARAM = MulticoreParam(6, 25, progressbar=TRUE)
gr <- countMAEReads(vcfs[1], rnas[1], BPPARAM=BPPARAM)

