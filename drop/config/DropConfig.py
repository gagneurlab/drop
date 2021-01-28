from .SampleAnnotation import SampleAnnotation
from .submodules import *
from .ExportCounts import ExportCounts
from drop import utils
from pathlib import Path
import wbuild


class DropConfig:
    CONFIG_KEYS = [
        # wbuild keys
        "projectTitle", "htmlOutputPath", "scriptsPath", "indexWithFolderName", "fileRegex", "readmePath",
        # global parameters
        "root", "sampleAnnotation", "geneAnnotation", "genomeAssembly", "exportCounts", "tools", "hpoFile",
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
        self.geneAnnotation = self.get("geneAnnotation")
        self.genomeAssembly = self.get("genomeAssembly")
        self.fastaFile = self.get("mae")["genome"] # TODO: move fasta outside of mae
        self.fastaDict = Path(self.fastaFile).with_suffix(".dict")
        self.sampleAnnotation = SampleAnnotation(self.get("sampleAnnotation"), self.root)

        # submodules
        self.AE = AE(
            config=self.get("aberrantExpression"),
            sampleAnnotation=self.sampleAnnotation,
            processedDataDir=self.processedDataDir,
            processedResultsDir=self.processedResultsDir
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
            processedResultsDir=self.processedResultsDir
        )

        # counts export
        self.exportCounts = ExportCounts(
            dict_=self.get("exportCounts"),
            outputRoot=self.processedResultsDir,
            sampleAnnotation=self.sampleAnnotation,
            geneAnnotations=self.getGeneAnnotations(),
            genomeAssembly=self.get("genomeAssembly"),
            aberrantExpression=self.AE,
            aberrantSplicing=self.AS
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

    def getGeneAnnotations(self):
        return self.geneAnnotation

    def getGeneVersions(self):
        return self.geneAnnotation.keys()

    def getGeneAnnotationFile(self, annotation):
        return self.geneAnnotation[annotation]

    def getFastaFile(self, str_=True):
        return utils.returnPath(self.fastaFile, str_)

    def getFastaDict(self, str_=True):
        return utils.returnPath(self.fastaDict, str_)

    def getBSGenomeName(self):
        assemblyID = self.get("genomeAssembly")

        if assemblyID == 'hg19':
            return "BSgenome.Hsapiens.UCSC.hg19"
        if assemblyID == 'hs37d5':
            return "BSgenome.Hsapiens.1000genomes.hs37d5"
        if assemblyID == 'hg38':
            return "BSgenome.Hsapiens.UCSC.hg38"
        if assemblyID == 'GRCh38':
            return "BSgenome.Hsapiens.NCBI.GRCh38"
        
        raise ValueError("Provided genome assembly not known: " + assemblyID)
 
    def getBSGenomeVersion(self):
        assemblyID = self.get("genomeAssembly")

        if assemblyID in ['hg19', 'hs37d5']:
            return 37
        if assemblyID in ['hg38', 'GRCh38']:
            return 38
        
        raise ValueError("Provided genome assembly not known: " + assemblyID)

    def getMafDbName(self):
        assemblyID = self.get("genomeAssembly")

        if assemblyID in ['hg19', 'hs37d5']:
            return "MafDb.gnomAD.r2.1.hs37d5"
        if assemblyID in ['hg38', 'GRCh38']:
            return "MafDb.gnomAD.r2.1.GRCh38"

        raise ValueError("Provided genome assembly not known: " + assemblyID)
 
 
