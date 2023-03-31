import numpy as np
import pandas as pd

from snakemake.io import expand

from drop import utils
from .Submodules import Submodule


class AS(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir, workDir):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir, workDir)
        self.CONFIG_KEYS = [
            "groups", "recount", "longRead", "filter", "minExpressionInOneSample", "minDeltaPsi",
            "quantileMinExpression", "quantileForFiltering", "implementation", "padjCutoff", 
            "deltaPsiCutoff", "maxTestedDimensionProportion", "genesToTest", "FRASER_version"
        ]
        self.name = "AberrantSplicing"
        # if self.run is false return without doing any config/sa checks for completeness
        if not self.run:
            return
        
        self.rnaIDs   = self.sampleAnnotation.subsetGroups(self.groups, assay="RNA")
        self.rnaExIDs = self.sampleAnnotation.subsetGroups(self.groups, assay="SPLICE_COUNT")
        for g in self.groups:
            if len(set(self.rnaIDs[g]) & set(self.rnaExIDs[g])) > 0:
                raise ValueError(f"{set(self.rnaIDs[g]) & set(self.extRnaIDs[g])} has both BAM and external count file \
                please fix in sample annotation table to only have either external count or BAM processing\n")

        all_ids = self.sampleAnnotation.subsetGroups(self.groups, assay=["RNA", "SPLICE_COUNT"])
        self.checkSubset(all_ids)

    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        setKey(dict_, None, "run", False)
        setKey(dict_, None, "groups", self.sampleAnnotation.getGroups(assay="RNA"))
        setKey(dict_, None, "recount", False)
        setKey(dict_, None, "longRead", False)
        setKey(dict_, None, "keepNonStandardChrs", False)
        setKey(dict_, None, "filter", True)
        setKey(dict_, None, "minExpressionInOneSample", 20)
        setKey(dict_, None, "quantileMinExpression", 10)
        setKey(dict_, None, "quantileForFiltering", 0.95)
        setKey(dict_, None, "minDeltaPsi", 0)
        setKey(dict_, None, "implementation", "PCA")
        setKey(dict_, None, "padjCutoff", 0.05)
        setKey(dict_, None, "deltaPsiCutoff", 0.05)
        setKey(dict_, None, "maxTestedDimensionProportion", 6)
        setKey(dict_, None, "genesToTest", None)
        setKey(dict_, None, "FRASER_version", "FRASER")
        return dict_

    def getSplitCountFiles(self, dataset):
        """
        Get all dummy count filenames for split counts
        :param dataset: DROP group name from wildcard
        :return: list of files
        """
        ids = self.sampleAnnotation.getIDsByGroup(dataset, assay="RNA")
        file_stump = self.processedDataDir / "aberrant_splicing" / "datasets" / "cache" / f"raw-local-{dataset}" / \
                     "sample_tmp" / "splitCounts"
        done_files = str(file_stump / "sample_{sample_id}.done")
        return expand(done_files, sample_id=ids)

    def getNonSplitCountFiles(self, dataset):
        """
        Get all dummy count filenames for non-split counts
        :param dataset: DROP group name from wildcard
        :return: list of files
        """
        ids = self.sampleAnnotation.getIDsByGroup(dataset, assay="RNA")
        file_stump = self.processedDataDir / "aberrant_splicing" / "datasets" / "cache" / f"raw-local-{dataset}" / \
                     "sample_tmp" / "nonSplitCounts"
        done_files = str(file_stump / "sample_{sample_id}.done")
        return expand(done_files, sample_id=ids)


    def getExternalCounts(self, group: str, fileType: str = "k_j_counts"):
        """
        Get externally provided splice count data dir based on the given group.
        If a file type is given the corresponding file within the folder is returned. 
        :param group: DROP group name from wildcard
        :param fileType: name of the file without extension which is to be returned
        :return: list of directories or files
        """

        # if sample annotation table does not contain SPLICE_COUNTS_DIR column. return no external counts
        if("SPLICE_COUNTS_DIR" not in self.sampleAnnotation.SAMPLE_ANNOTATION_COLUMNS):
            return []

        ids = self.sampleAnnotation.getIDsByGroup(group, assay="SPLICE_COUNT")
        extCountFiles = self.sampleAnnotation.getImportCountFiles(annotation=None, group=group, 
                file_type="SPLICE_COUNTS_DIR", asSet=False)
        if fileType is not None:
            extCountFiles = np.asarray(extCountFiles)[pd.isna(extCountFiles) == False].tolist()
            extCountFiles = [x + "/" + fileType + ".tsv.gz" for x in extCountFiles]
        return extCountFiles

    def getPsiTypeAssay(self):    
        """
        Dependent on the FRASER version, get the psiType of the h5 assays that 
        will be produced by the aberrant-splicing pipeline.
        :return: psiType
        """
        
        fraser_version = self.get("FRASER_version")
        if(fraser_version == "FRASER2"):
            return "jaccard"
        return "theta"
