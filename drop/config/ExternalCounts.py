from snakemake.io import expand

class ExternalCounts:

    def __init__(self, dropConfig):
        self.cfg = dropConfig
        self.sa = dropConfig.sampleAnnotation
        self.processedResults = self.cfg.getProcessedResultsDir(str_=False)
        self.COUNT_TYPE_MAP = {
            "geneCounts": "aberrantExpression",
            "splitCounts": "aberrantSplicing",
            "spliceSiteOverlapCounts": "aberrantSplicing"
        }

    def getExportGroups(self, modules=None):
        if modules is None:
            modules = ["aberrantExpression", "aberrantSplicing"]
        groups = []  # get all active groups
        for module in modules:
            groups.extend(self.cfg.get(module)["groups"])

        exclude = self.cfg.get("exportCounts")["excludeGroups"]
        return set(groups) - set(exclude)

    def getExportCountFiles(self, prefix):
        if prefix not in self.COUNT_TYPE_MAP.keys():
            raise ValueError(f"{prefix} not a valid file type for exported counts")

        datasets = self.getExportGroups([self.COUNT_TYPE_MAP[prefix]])
        annotations = self.cfg.get("exportCounts")["geneAnnotations"]
        genomeAssembly = self.cfg.get("genomeAssembly")

        pattern = self.processedResults / "exported_counts" / "{dataset}--{genomeAssembly}--{annotation}" / f"{prefix}.tsv.gz"
        return expand(str(pattern), annotation=annotations, dataset=datasets, genomeAssembly=genomeAssembly)

    def getImportCountFiles(self, annotation, group, file_type="GENE_COUNTS_FILE",
                            annotation_key="ANNOTATION", group_key="DROP_GROUP"):
        """
        :param annotation: annotation name as specified in config and ANNOTATION column
        :param group: a group of the DROP_GROUP column
        :return: set of unique external count file names
        """
        subset = self.sa.subsetSampleAnnotation(annotation_key, annotation)
        subset = self.sa.subsetSampleAnnotation(group_key, group, subset)
        return set(subset[file_type].tolist())

