from .SampleAnnotation import SampleAnnotation
from .Submodules import AE, AS, MAE
from drop import utils
import wbuild
from pathlib import Path
from snakemake.logging import logger
from snakemake.io import expand
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')

class DropConfig:
    
    FILE_KEYS = ["htmlOutputPath", "root", "geneAnnotation", "sampleAnnotation", "mae"]
    CONFIG_KEYS = [] # TODO: list all keys
    
    def __init__(self, config):
        
        config = self.checkConfig(config)
        self.config = self.setDefaults(config)
        
        self.root = Path(self.get("root"))
        self.processedDataDir = self.root / "processed_data"
        self.processedResultsDir = self.root / "processed_results"
        self.htmlOutputPath = Path(self.get("htmlOutputPath"))
        utils.createIfMissing(self.root)
        utils.createIfMissing(self.processedDataDir)
        utils.createIfMissing(self.processedResultsDir)
        
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
        return utils.returnPath(self.processedResultsDir)
    
    def getHtmlOutputPath(self, str_=True):
        return utils.returnPath(self.htmlOutputPath)
    
    def getHtmlFromScript(self, path):
        stump = self.htmlOutputPath / utils.getRuleFromPath(path, prefix=True)
        return str(stump) + ".html"
    
    def get(self, key):
        if key not in self.config:
            if key not in self.CONFIG_KEYS:
                raise KeyError(f"{key} not defined for Drop config")
            raise KeyError(f"{key} missing in config dictionary")
        return self.config[key]
    
    def getGeneAnnotations(self):
        return self.geneAnnotation
        
    def getGeneVersions(self):
        return self.geneAnnotation.keys()
    
    def getGeneAnnotationFile(self, annotation):
        return self.geneAnnotation[annotation]
    
    def getAE(self):
        return self.AE
    
    def getAS(self):
        return self.AS
    
    def getMAE(self):
        return self.MAE
    
    # TODO: wrap in separate class
    def getExportGroups(self, modules=None):
        if modules is None:
            modules = ["aberrantExpression", "aberrantSplicing"]
        groups = [] # get all active groups
        for module in modules:
            groups.extend(self.get(module)["groups"])
        exclude = self.get("exportCounts")["excludeGroups"]
        return set(groups) - set(exclude)
        
    def getExportCountFiles(self, prefix):

        count_type_map = {"geneCounts":"aberrantExpression",
                          "splitCounts":"aberrantSplicing",
                          "spliceSiteOverlapCounts":"aberrantSplicing"}
        if prefix not in count_type_map.keys():
            raise ValueError(f"{prefix} not a valid file type for exported counts")

        datasets = self.getExportGroups([count_type_map[prefix]])
        annotations = self.get("exportCounts")["geneAnnotations"]
        genomeAssembly = self.get("genomeAssembly")

        pattern = self.getProcessedResultsDir() + f"/exported_counts/{{dataset}}--{{genomeAssembly}}--{{annotation}}/{prefix}.tsv.gz"
        return expand(pattern, annotation=annotations, dataset=datasets, genomeAssembly=genomeAssembly)
    

