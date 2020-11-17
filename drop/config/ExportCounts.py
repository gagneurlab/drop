from snakemake.io import expand
from drop import utils

class ExportCounts:

    COUNT_TYPE_MAP = {
        "geneCounts": "aberrantExpression",
        "splicingCounts": "aberrantSplicing",
    }

    def __init__(self, dict_, outputRoot, sampleAnnotation, geneAnnotations, genomeAssembly,
                 aberrantExpression, aberrantSplicing):
        """
        :param dict_: config dictionary for count export
        :param sampleAnnotation: parsed sample annotation
        :param geneAnnotations: list of gene annotation names
        :param aberrantExpression: AberrantExpression object
        :param aberrantSplicing: AberrantSplicing object
        """
        self.CONFIG_KEYS = ["geneAnnotations", "excludeGroups"]
        self.config_dict = self.setDefaults(dict_, geneAnnotations)
        self.outputRoot = outputRoot
        self.sa = sampleAnnotation
        self.genomeAssembly = genomeAssembly
        self.modules = {
            "aberrantExpression": aberrantExpression,
            "aberrantSplicing": aberrantSplicing
        }
        self.pattern = self.outputRoot / "exported_counts" / "{dataset}--{genomeAssembly}--{annotation}"

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
            str_=True
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
        return export_groups

    def getExportCountFiles(self, prefix, expandPattern=None, **kwargs):
        """
        Determine export count files.
        :param prefix: name of file
        :return: list of files to
        """
        if prefix not in self.COUNT_TYPE_MAP.keys():
            raise ValueError(f"{prefix} not a valid file type for exported counts")

        datasets = self.getExportGroups([self.COUNT_TYPE_MAP[prefix]])
        if expandPattern is None:
            file_pattern = str(self.pattern / f"{prefix}.tsv.gz")
        else:
            file_pattern = str(self.pattern / f"{expandPattern}.tsv.gz")
        count_files = expand(file_pattern, annotation=self.get("geneAnnotations"),
                             dataset=datasets, genomeAssembly=self.genomeAssembly, **kwargs)
        return count_files

