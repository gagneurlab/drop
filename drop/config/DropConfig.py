from .SampleAnnotation import SampleAnnotation
from .Submodules import AE, AS, MAE
from .ExternalCounts import ExternalCounts
from drop import utils
from pathlib import Path

class DropConfig:
    
    FILE_KEYS = ["htmlOutputPath", "root", "geneAnnotation", "sampleAnnotation", "mae"]
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
        config = self.checkConfig(wbuildConfig.getConfig())
        self.config = self.setDefaults(config)
        
        self.root = Path(self.get("root"))
        self.processedDataDir = self.root / "processed_data"
        self.processedResultsDir = self.root / "processed_results"
        utils.createIfMissing(self.root)
        utils.createIfMissing(self.processedDataDir)
        utils.createIfMissing(self.processedResultsDir)

        self.htmlOutputPath = Path(self.get("htmlOutputPath"))
        self.readmePath = Path(self.get("readmePath"))

        self.sampleAnnotation = SampleAnnotation(self.get("sampleAnnotation"), self.root)
        self.geneAnnotation = self.get("geneAnnotation")
        self.genomeAssembly = self.get("genomeAssembly")
        
        # setup submodules
        cfg = self.config
        sa = self.sampleAnnotation
        pd = self.processedDataDir
        pr = self.processedResultsDir
        self.AE = AE(cfg, sa, pd, pr)
        self.AS = AS(cfg, sa, pd, pr)
        self.MAE = MAE(cfg, sa, pd, pr)

        # external counts settings
        self.externalCounts = ExternalCounts(self)
        
        # legacy
        utils.setKey(self.config, None, "aberrantExpression", self.AE.dict_)
        utils.setKey(self.config, None, "aberrantSplicing", self.AS.dict_)
        utils.setKey(self.config, None, "mae", self.MAE.dict_)
        
    def checkConfig(self, config):
        # TODO: check if files exists too!
        self.checkKeys(config, keys=self.FILE_KEYS)
        self.checkKeys(config["geneAnnotation"], keys=None)
        self.checkKeys(config["mae"], keys=["genome", "qcVcf"])
        return config
    
    def checkKeys(self, dict_, keys):
        keys = dict_.keys() if keys is None else keys
        for key in keys:
            if key not in dict_.keys():
                raise KeyError(f"{key} is mandatory but missing")
            # get real path
            if isinstance(dict_[key], str):
                filename = dict_[key]
                dict_[key] = str(Path(filename).expanduser())
    
    def setDefaults(self, config, method=None):
        """
        set defaults for config keys
        """
        config["indexWithFolderName"] = True
        config["fileRegex"] = ".*\.R"
        config["wBuildPath"] =  utils.getWBuildPath()
        
        setKey = utils.setKey
        setKey(config, None, "projectTitle", "DROP: Detection of RNA Outlier Pipeline")
        
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