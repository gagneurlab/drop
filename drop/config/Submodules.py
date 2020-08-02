from pathlib import Path
from snakemake.logging import logger
from snakemake.io import expand
from drop import utils

class Submodule:

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir):
        self.CONFIG_KEYS = []
        self.name = "Submodule"
        self.processedDataDir = processedDataDir
        self.processedResultsDir = processedResultsDir
        self.sa = sampleAnnotation
        self.dict_ = self.setDefaultKeys(config)
        self.groups = self.dict_["groups"]
    
    def setDefaultKeys(self, dict_):
        dict_ = {} if dict_ is None else dict_
        return dict_

    def get(self, key):
        if key not in self.CONFIG_KEYS:
            raise KeyError(f"{key} not defined for {self.name} config")
        return self.dict_[key]

    def getWorkdir(self, str_=True):
        return utils.returnPath(Path("Scripts") / self.name / "pipeline", str_)

class AE(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir, externalCounts):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir)
        self.CONFIG_KEYS = [
            "groups", "fpkmCutoff", "implementation", "padjCutoff", "zScoreCutoff",
            "maxTestedDimensionProportion"
        ]
        self.name = "AberrantExpression"
        self.extCounts = externalCounts
        self.rnaIDs = self.sa.subsetGroups(self.groups, assay="RNA")
        self.extRnaIDs = self.sa.subsetGroups(self.groups, assay="GENE_COUNTS", warn=0, error=0)
    
    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        setKey(dict_, None, "groups", self.sa.getGroups(assay="RNA"))
        setKey(dict_, None, "fpkmCutoff", 1)
        setKey(dict_, None, "implementation", "autoencoder")
        setKey(dict_, None, "padjCutoff", .05)
        setKey(dict_, None, "zScoreCutoff", 0)
        setKey(dict_, None, "maxTestedDimensionProportion", 3)
        return dict_

    def getCountFiles(self, annotation, group):
        """
        Get all count files from DROP (counted from BAM file) and external count matrices
        :param annotation: annotation name from wildcard
        :param group: DROP group name from wildcard
        :return: list of files
        """
        bam_IDs = self.sa.getIDsByGroup(group, assay="RNA")
        file_stump = self.processedDataDir / "aberrant_expression" / annotation / "counts" / "{sampleID}.Rds"
        count_files = expand(str(file_stump), sampleID=bam_IDs)
        extCountFiles = self.extCounts.getImportCountFiles(annotation, group, file_type="GENE_COUNTS_FILE")
        count_files.extend(extCountFiles)
        return count_files

    def getCountParams(self, rnaID):
        sa_row = self.sa.getRow("RNA_ID", rnaID)
        count_params = sa_row[["STRAND", "COUNT_MODE", "PAIRED_END", "COUNT_OVERLAPS"]]
        return count_params.iloc[0].to_dict()


class AS(Submodule):
    
    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir, externalCounts):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir)
        self.CONFIG_KEYS = [
            "groups", "recount", "longRead", "filter", "minExpressionInOneSample", "minDeltaPsi",
            "implementation", "padjCutoff", "zScoreCutoff", "deltaPsiCutoff", "maxTestedDimensionProportion"
        ]
        self.name = "AberrantSplicing"
        self.extCounts = externalCounts
        self.rnaIDs = self.sa.subsetGroups(self.groups, assay="RNA")
        
    
    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        setKey(dict_, None, "groups", self.sa.getGroups(assay="RNA"))
        setKey(dict_, None, "recount", False)
        setKey(dict_, None, "longRead", False)
        setKey(dict_, None, "filter", True)
        setKey(dict_, None, "minExpressionInOneSample", 20)
        setKey(dict_, None, "minDeltaPsi", 0)
        setKey(dict_, None, "implementation", "PCA")
        setKey(dict_, None, "padjCutoff", 0.05)
        setKey(dict_, None, "zScoreCutoff", 0.05)
        setKey(dict_, None, "deltaPsiCutoff", 0.05)
        setKey(dict_, None, "maxTestedDimensionProportion", 6)
        return dict_
    
class MAE(Submodule):
    
    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir)
        self.CONFIG_KEYS = [
            "groups", "genome", "qcVcf", "qcGroups", "gatkIgnoreHeaderCheck", "padjCutoff",
            "allelicRatioCutoff", "maxAF", "gnomAD"
        ]
        self.name = "MonoallelicExpression"
        self.qcGroups = self.dict_["qcGroups"]
        self.maeIDs = self.createMaeIDS(id_sep='--')
    
    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        dict_ = utils.checkKeys(dict_, keys=["genome", "qcVcf"], check_files=True)
        groups = setKey(dict_, None, "groups", self.sa.getGroups(assay="DNA"))
        setKey(dict_, None, "qcGroups", groups)
        setKey(dict_, None, "gatkIgnoreHeaderCheck", True)
        setKey(dict_, None, "padjCutoff", .05)
        setKey(dict_, None, "allelicRatioCutoff", 0.8)
        setKey(dict_, None, "maxAF", .001)
        setKey(dict_, None, "gnomAD", False)
        return dict_
    
    def createMaeIDS(self, id_sep='--'):
        """
        Create MAE IDs from sample annotation
        """
        grouped_rna_ids = self.sa.subsetGroups(self.groups, assay="RNA", warn=1, error=1)
        id_map = self.sa.idMapping
        mae_ids = {}
        for gr, rna_ids in grouped_rna_ids.items():
            subset = id_map[id_map["RNA_ID"].isin(rna_ids)]
            dna_rna_pairs = zip(subset["DNA_ID"], subset["RNA_ID"])
            mae_ids[gr] = list(map(id_sep.join, dna_rna_pairs))
        return mae_ids
    
    def getMaeByGroup(self, group):
        if not isinstance(group, str):
            group = list(group)[0]
        return self.maeIDs[group]

    def getMaeAll(self):
        """
        Get a list of all MAE IDs from the groups specified in the config.
        Useful for collecting all MAE IDs ungrouped.
        """
        all_ids = []
        groups = [self.groups] if isinstance(self.groups, str) else self.groups
        for group in groups:
            all_ids.extend(self.getMaeByGroup(group))
        return all_ids
