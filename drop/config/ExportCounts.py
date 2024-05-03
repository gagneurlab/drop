from snakemake.io import expand
from drop import utils
from snakemake.logging import logger


class ExportCounts:
    COUNT_TYPE_MAP = {
        "geneCounts": "aberrantExpression",
        "splicingCounts": "aberrantSplicing",
    }

    def __init__(
            self,
            dict_,
            outputRoot,
            sampleAnnotation,
            genome,
            aberrantExpression,
            aberrantSplicing
    ):
        """
        :param dict_: config dictionary for count export
        :param sampleAnnotation: parsed sample annotation
        :param geneAnnotations: list of gene annotation names
        :param aberrantExpression: AberrantExpression object
        :param aberrantSplicing: AberrantSplicing object
        """
        self.CONFIG_KEYS = ["geneAnnotations", "excludeGroups"]
        self.config_dict = self.setDefaults(dict_, genome.annotation)
        self.outputRoot = outputRoot / "exported_counts"
        self.sampleAnnotation = sampleAnnotation
        self.genomeAssembly = genome.assembly
        self.geneAnnotations = self.get("geneAnnotations")
        self.modules = {
            "aberrantExpression": aberrantExpression,
            "aberrantSplicing": aberrantSplicing
        }

        self.pattern = self.outputRoot / "{dataset}--{genomeAssembly}--{annotation}"

        self.checkNonExternalGeneAnnotation()

    def setDefaults(self, config_dict, gene_annotations):
        utils.setKey(config_dict, None, "geneAnnotations", gene_annotations)
        utils.setKey(config_dict, None, "excludeGroups", list())

        # check consistency of gene annotations
        anno_incomp = set(config_dict["geneAnnotations"]) - set(gene_annotations)
        if len(anno_incomp) > 0:
            message = f"{anno_incomp} are not valid annotation version in 'geneAnnotation'"
            message += "but required in 'exportCounts'.\n Please make sure they match."
            raise ValueError(message)

        return config_dict

    def get(self, key):
        if key not in self.CONFIG_KEYS:
            raise KeyError(f"{key} not defined for count export")
        return self.config_dict[key]

    def getFilePattern(self, str_=True, expandStr=False):
        pattern = self.pattern
        if expandStr:
            pattern = pattern.__str__().replace("{", "{{").replace("}", "}}")
        return utils.returnPath(pattern, str_=str_)

    def getExportGroups(self, modules=None):
        """
        Determine from which DROP groups counts should be exported
        :param modules: 'aberrantExpression' for gene counts, 'aberrantSplicing' for splicing counts export
        :return: DROP groups from which to export counts
        """
        if modules is None:
            modules = self.modules.keys()
        elif isinstance(modules, str):
            modules = [modules]
        groups = []  # get all active groups
        for module in modules:
            groups.extend(self.modules[module].groups)
        export_groups = set(groups) - set(self.get("excludeGroups"))

        return sorted(list(export_groups))

    def getFiles(self, filename, datasets=None, **kwargs):
        """
        Determine files for export count groups.
        :param filename: name of file
        :return: list of export files
        """
        if datasets is None:
            datasets = self.getExportGroups()
        file_pattern = str(self.pattern / f"{filename}")
        return expand(
            file_pattern,
            dataset=datasets,
            annotation=self.geneAnnotations,
            genomeAssembly=self.genomeAssembly,
            **kwargs
        )

    def getExportCountFiles(self, count_type, suffix="tsv.gz", expandPattern=None, **kwargs):
        """
        Determine export count files.
        :param count_type: count type for mapping the submodule
        :param suffix: file type suffix (without dot)
        :return: list of export count files
        """
        if count_type not in self.COUNT_TYPE_MAP.keys():
            raise ValueError(f"'{count_type}' not a valid file type for exported counts")
        datasets = self.getExportGroups([self.COUNT_TYPE_MAP[count_type]])
        expandPattern = count_type if expandPattern is None else expandPattern
        return self.getFiles(f"{expandPattern}.{suffix}", datasets, **kwargs)

    def checkNonExternalGeneAnnotation(self):
        excluded_groups = self.config_dict['excludeGroups']
        print(excluded_groups)

        non_excluded_samples = self.sampleAnnotation.annotationTable[self.sampleAnnotation.annotationTable['DROP_GROUP'].isin(excluded_groups) == False]
        print(non_excluded_samples)
        if sum(non_excluded_samples['GENE_ANNOTATION'].isna() == False) > 0:
            logger.info("WARNING: Found %d samples that had `GENE_ANNOTATION` provided in sample annotation table but are not external counts. The provided `GENE_ANNOTATIONs` are ignored.\n" % (sum(non_excluded_samples['GENE_ANNOTATION'].isna() == False)))
            self.sampleAnnotation.annotationTable.loc[self.sampleAnnotation.annotationTable['DROP_GROUP'].isin(excluded_groups) == False, "GENE_ANNOTATION"] = ""
