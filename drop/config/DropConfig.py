from .SampleAnnotation import SampleAnnotation
from .Submodules import AE, AS, MAE
from .ExternalCounts import ExternalCounts
from drop import utils
from pathlib import Path

class DropConfig:

    CONFIG_KEYS = [
        # wbuild keys
        "projectTitle", "htmlOutputPath", "scriptsPath", "indexWithFolderName", "fileRegex", "readmePath",
        # global parameters
        "root", "sampleAnnotation", "geneAnnotation", "genomeAssembly", "scanBamParam", "exportCounts", "tools",
        # modules
        "aberrantExpression", "aberrantSplicing", "mae"
    ]
    
    def __init__(self, wbuildConfig):
        """
        Parse wbuild/snakemake config object for DROP-specific content

        :param wbuildConfig: wBuild config object
        """

        self.wBuildConfig = wbuildConfig
        self.config = self.setDefaults(wbuildConfig.getConfig())
        
        self.root = Path(self.get("root"))
        self.processedDataDir = self.root / "processed_data"
        self.processedResultsDir = self.root / "processed_results"
        utils.createDir(self.root)
        utils.createDir(self.processedDataDir)
        utils.createDir(self.processedResultsDir)

        self.htmlOutputPath = Path(self.get("htmlOutputPath"))
        self.readmePath = Path(self.get("readmePath"))

        self.geneAnnotation = self.get("geneAnnotation")
        self.genomeAssembly = self.get("genomeAssembly")
        self.sampleAnnotation = SampleAnnotation(self.get("sampleAnnotation"), self.root)
        self.externalCounts = ExternalCounts(self)

        # setup submodules
        cfg = self.config
        sa = self.sampleAnnotation
        pd = self.processedDataDir
        pr = self.processedResultsDir
        ec = self.externalCounts
        self.AE = AE(cfg["aberrantExpression"], sa, pd, pr, ec)
        self.AS = AS(cfg["aberrantSplicing"], sa, pd, pr, ec)
        self.MAE = MAE(cfg["mae"], sa, pd, pr)

        # legacy
        utils.setKey(self.config, None, "aberrantExpression", self.AE.dict_)
        utils.setKey(self.config, None, "aberrantSplicing", self.AS.dict_)
        utils.setKey(self.config, None, "mae", self.MAE.dict_)


    def setDefaults(self, config):
        """
        Check mandatory keys and set defaults for any missing keys
        :param config: config dictionary
        :return: config dictionary with defaults
        """
        # check mandatory keys
        config = utils.checkKeys(config, keys=["htmlOutputPath", "root", "sampleAnnotation"], check_files=True)
        config["geneAnnotation"] = utils.checkKeys(config["geneAnnotation"], keys=None, check_files=True)

        config["indexWithFolderName"] = True
        config["fileRegex"] = ".*\.R"
        config["wBuildPath"] = utils.getWBuildPath()
        
        setKey = utils.setKey
        setKey(config, None, "genomeAssembly", "hg19")
        setKey(config, None, "scanBamParam", "null")
        
        # export settings
        setKey(config, None, "exportCounts", dict())
        gene_annotations = list(config["geneAnnotation"].keys())
        setKey(config, ["exportCounts"], "geneAnnotation", gene_annotations)
        setKey(config, ["exportCounts"], "excludeGroups", list())
        
        # check consistency of gene annotations
        anno_incomp = set(config["exportCounts"]["geneAnnotations"]) - set(gene_annotations)
        if len(anno_incomp) > 0:
            message = f"{anno_incomp} are not valid annotation version in 'geneAnnotation'"
            message += "but required in 'exportCounts'.\n Please make sure they match."
            raise ValueError(message)
        
        # set submodule dictionaries
        setKey(config, None, "aberrantExpression", dict())
        setKey(config, None, "aberrantSplicing", dict())
        setKey(config, None, "mae", dict())
        
        # commandline tools
        setKey(config, None, "tools", dict())
        setKey(config, ["tools"], "samtoolsCmd", "samtools")
        setKey(config, ["tools"], "bcftoolsCmd", "bcftools")
        setKey(config, ["tools"], "gatkCmd", "gatk")
        
        return config
    
    def getRoot(self, str_=True):
        return utils.returnPath(self.root, str_=str_)
    
    def getProcessedDataDir(self, str_=True):
        return utils.returnPath(self.processedDataDir, str_=str_)
    
    def getProcessedResultsDir(self, str_=True):
        return utils.returnPath(self.processedResultsDir, str_=str_)
    
    def getHtmlOutputPath(self, str_=True):
        return utils.returnPath(self.htmlOutputPath, str_=str_)

    #def getReadmePath(self, str_=True):
    #    readme_name = Path(self.readmePath).name
    #    readme_name = readme_name.replace(".md", ".html")
    #    return utils.returnPath(self.htmlOutputPath / readme_name, str_=str_)
    
    def getHtmlFromScript(self, path):
        stump = self.htmlOutputPath / utils.getRuleFromPath(path, prefix=True)
        return str(stump) + ".html"
    
    def get(self, key):
        if key not in self.CONFIG_KEYS:
            raise KeyError(f"{key} not defined for Drop config")
        return self.wBuildConfig.get(key)

    def getGeneAnnotations(self):
        return self.geneAnnotation
        
    def getGeneVersions(self):
        return self.geneAnnotation.keys()
    
    def getGeneAnnotationFile(self, annotation):
        return self.geneAnnotation[annotation]
