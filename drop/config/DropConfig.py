from .SampleAnnotation import SampleAnnotation
from .Genome import Genome
from .SampleParams import SampleParams
from .submodules import *
from .ExportCounts import ExportCounts
from drop import utils
from pathlib import Path
import wbuild
from snakemake.logging import logger


class DropConfig:
    CONFIG_KEYS = [
        # wbuild keys
        "projectTitle", "htmlOutputPath", "scriptsPath", "indexWithFolderName", "fileRegex", "readmePath",
        # global parameters
        "root", "sampleAnnotation", "geneAnnotation", "genomeAssembly", "exportCounts", "tools", "hpoFile","genome",
        # modules
        "aberrantExpression", "aberrantSplicing", "mae"
    ]

    def __init__(self, wbuildConfig):
        """
        Parse wbuild/snakemake config object for DROP-specific content
        :param wbuildConfig: wBuild config object
        """

        self.wBuildConfig = wbuildConfig
        self.config_dict = self.setDefaults(wbuildConfig.getConfig())

        self.root = Path(self.get("root"))
        self.processedDataDir = self.root / "processed_data"
        self.processedResultsDir = self.root / "processed_results"
        utils.createDir(self.root)
        utils.createDir(self.processedDataDir)
        utils.createDir(self.processedResultsDir)

        self.htmlOutputPath = Path(self.get("htmlOutputPath"))
        self.readmePath = Path(self.get("readmePath"))

        # annotations
        self.genome = Genome(
            annotation=self.get("geneAnnotation"),
            assembly=self.get("genomeAssembly"),
            reference=self.get("genome")
        )

        self.sampleAnnotation = SampleAnnotation(
            file=self.get("sampleAnnotation"),
            root=self.root,
            genome=self.genome
        )

        # submodules
        self.AE = AE(
            config=self.get("aberrantExpression"),
            sampleAnnotation=self.sampleAnnotation,
            processedDataDir=self.processedDataDir,
            processedResultsDir=self.processedResultsDir,
        )

        self.AS = AS(
            config=self.get("aberrantSplicing"),
            sampleAnnotation=self.sampleAnnotation,
            processedDataDir=self.processedDataDir,
            processedResultsDir=self.processedResultsDir
        )
        self.MAE = MAE(
            config=self.get("mae"),
            sampleAnnotation=self.sampleAnnotation,
            processedDataDir=self.processedDataDir,
            processedResultsDir=self.processedResultsDir,
            genome=self.genome
        )

        # counts export
        self.exportCounts = ExportCounts(
            dict_=self.get("exportCounts"),
            outputRoot=self.processedResultsDir,
            sampleAnnotation=self.sampleAnnotation,
            genome=self.genome,
            aberrantExpression=self.AE,
            aberrantSplicing=self.AS
        )

        # write sample params for each module AS not currently supported
        sampleParams = SampleParams(
            self.AE, 
            self.MAE, 
            self.get("geneAnnotation"),
            self.processedDataDir, 
            self.sampleAnnotation
        )


        # legacy
        utils.setKey(self.config_dict, None, "aberrantExpression", self.AE.dict_)
        utils.setKey(self.config_dict, None, "aberrantSplicing", self.AS.dict_)
        utils.setKey(self.config_dict, None, "mae", self.MAE.dict_)

    def setDefaults(self, config_dict):
        """
        Check mandatory keys and set defaults for any missing keys
        :param config_dict: config dictionary
        :return: config dictionary with defaults
        """
        # check mandatory keys
        config_dict = utils.checkKeys(config_dict, keys=["htmlOutputPath", "root", "sampleAnnotation"],
                                      check_files=True)
        config_dict["geneAnnotation"] = utils.checkKeys(config_dict["geneAnnotation"], keys=None, check_files=True)

        config_dict["wBuildPath"] = utils.getWBuildPath()

        setKey = utils.setKey
        setKey(config_dict, None, "fileRegex", r".*\.(R|md)")
        setKey(config_dict, None, "genomeAssembly", "hg19")
        hpo_url = 'https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/hpo_genes.tsv.gz'
        setKey(config_dict, None, "hpoFile", hpo_url)

        # set submodule dictionaries
        setKey(config_dict, None, "aberrantExpression", dict())
        setKey(config_dict, None, "aberrantSplicing", dict())
        setKey(config_dict, None, "mae", dict())
        setKey(config_dict, None, "exportCounts", dict())

        # Legacy check: If mae still defines genome print warning, otherwise use the
        # globally defined genome
        try:
            genome_files = self.get("mae")["genome"]
            logger.info(
                "WARNING: Using the mae defined genome instead of the globally defined one.\n"
                "This will be deprecated in the future to allow for reference genomes to"
                " be defined in the sample annotation table. Please update your config "
                "and sample annotation table\n"
            )
        except KeyError:
            genome_files = self.get("genome")
        setKey(config_dict, None, "genome", genome_files)


        # commandline tools
        setKey(config_dict, None, "tools", dict())
        setKey(config_dict, ["tools"], "samtoolsCmd", "samtools")
        setKey(config_dict, ["tools"], "bcftoolsCmd", "bcftools")
        setKey(config_dict, ["tools"], "gatkCmd", "gatk")

        return config_dict

    def getRoot(self, str_=True):
        return utils.returnPath(self.root, str_=str_)

    def getProcessedDataDir(self, str_=True):
        return utils.returnPath(self.processedDataDir, str_=str_)

    def getProcessedResultsDir(self, str_=True):
        return utils.returnPath(self.processedResultsDir, str_=str_)

    def getHtmlOutputPath(self, str_=True):
        return utils.returnPath(self.htmlOutputPath, str_=str_)

    def getHtmlFromScript(self, path, str_=True):
        path = Path(path).with_suffix(".html")
        file_name = wbuild.utils.pathsepsToUnderscore(str(path), dotsToUnderscore=False)
        html_output_path = self.getHtmlOutputPath(str_=False)
        return utils.returnPath(html_output_path / file_name, str_=str_)

    def get(self, key):
        if key not in self.CONFIG_KEYS:
            raise KeyError(f"'{key}' not defined for DROP config")
        return self.wBuildConfig.get(key)

    def getTool(self, tool):
        try:
            toolCmd = self.get("tools")[tool]
        except KeyError:
            raise KeyError(f"'{toolCmd}' not a defined tool for DROP config")
        return toolCmd
