import itertools

cfg = snakemake.params.dropConfig
exportCounts = cfg.exportCounts
sampleAnnotation = snakemake.params.sampleAnnotation

# adapt sample annotation data frame
sa = sampleAnnotation.annotationTable.copy()
sa["STRAND_SPECIFIC"] = sa["STRAND"] != "no"
if "ICD_10" in sa:
    sa["ICD_10"] = sa["ICD_10"].str.upper()

sa_cols = frozenset(['RNA_ID', 'INDIVIDUAL_ID', 'TISSUE', 'SEX', 'AFFECTED', 'ICD_10', 'PAIRED_END',
                     'STRAND_SPECIFIC']).intersection(sa.columns)
sa_cols = list(sa_cols)


# Getters for DESCRIPTION file
def get_tissue_info(sa):
    if "TISSUE" not in sa or sa["TISSUE"].isnull().all():
            return "No tissue information provided"
    uniq_tissues = sa["TISSUE"].unique()
    if len(uniq_tissues) > 1:
        raise ValueError("Either missing value or more than 1 tissue in dataset")
    return "".join(uniq_tissues)


def get_disease_info(sa):
    if "ICD_10" in sa:
        # unite(data.table(table(sa_sub$ICD_10)), col = 'aux', 'V1', 'N', sep = ': ')$aux, collapse = ', ' )
        table = sa.ICD_10.value_counts()
        return "\n" + "\n".join([f"\t{d}: {c}" for d, c in zip(table.index, table)])
    return "No disease information"


def get_strand(sa):
    if sa.STRAND_SPECIFIC.nunique() > 1:
        return "All samples should be either strand- or non-strand-specific!"
    else:
        return sa.STRAND_SPECIFIC.unique()[0]


def get_pairedEnd(sa):
    if sa.PAIRED_END.nunique() > 1:
        return 'All samples should be either single end or paired end!',
    return sa.PAIRED_END.unique()[0]


desc = """Title: # Add a title
Number of samples: {}
Tissue: {}
Organism: {}
Genome assembly: {}
Gene annotation: {}
Disease (ICD-10: N): {}
Strand specific: {}
Paired end: {}
Cite as: RNA-Seq count tables were taken from # add your citation(s)
Dataset contact: # Use format Name Last_Name, <email address>
Comments: # add any comments, if needed, otherwise remove
"""

group_with_anno = itertools.product(exportCounts.getExportGroups(), exportCounts.geneAnnotations)
sa_files = snakemake.output.sampleAnnotations
desc_files = snakemake.output.descriptions

for ga, sa_file, desc_file in zip(group_with_anno, sa_files, desc_files):
    group, annotation = ga # unpack

    # subset by group
    rna_ids = sampleAnnotation.rnaIDs[group]
    sa_sub = sa.loc[sa["RNA_ID"].isin(rna_ids)]

    # save sample annotation subset
    sa_sub.to_csv(sa_file, sep='\t', index=False, columns=sa_cols)

    # save DESCRIPTION file
    with open(desc_file, "w") as f:
        desc_output = desc.format(
            len(sa_sub),  # number samples
            get_tissue_info(sa_sub),  # tissue
            "Homo sapiens",  # organism
            exportCounts.genomeAssembly,
            annotation,
            get_disease_info(sa_sub),  # disease
            get_strand(sa_sub),  # strand specific
            get_pairedEnd(sa_sub)  # paired end
        )
        f.write(desc_output)
