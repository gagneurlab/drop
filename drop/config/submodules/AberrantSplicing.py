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
            "implementation", "padjCutoff", "zScoreCutoff", "deltaPsiCutoff", "maxTestedDimensionProportion"
        ]
        self.name = "AberrantSplicing"
        # if self.run is false return without doing any config/sa checks for completeness
        if not self.run:
            return
        self.rnaIDs = self.sampleAnnotation.subsetGroups(self.groups, assay="RNA")
        self.checkSubset(self.rnaIDs)

    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        setKey(dict_, None, "groups", self.sampleAnnotation.getGroups(assay="RNA"))
        setKey(dict_, None, "recount", False)
        setKey(dict_, None, "longRead", False)
        setKey(dict_, None, "keepNonStandardChrs", False)
        setKey(dict_, None, "filter", True)
        setKey(dict_, None, "minExpressionInOneSample", 20)
        setKey(dict_, None, "minDeltaPsi", 0)
        setKey(dict_, None, "implementation", "PCA")
        setKey(dict_, None, "padjCutoff", 0.05)
        setKey(dict_, None, "zScoreCutoff", 0.05)
        setKey(dict_, None, "deltaPsiCutoff", 0.05)
        setKey(dict_, None, "maxTestedDimensionProportion", 6)
        return dict_

    def getSplitCountFiles(self, dataset):
        """
        Get all dummy count filenames for split counts
        :param dataset: DROP group name from wildcard
        :return: list of files
        """
        ids = self.sampleAnnotation.getIDsByGroup(dataset, assay="RNA")
        file_stump = self.processedDataDir / "aberrant_splicing" / "datasets" / "cache"
        count_file = str(file_stump / "splitCounts-{sample_id}.RDS")
        return expand(count_file, sample_id=ids)


    def getNonSplitCountFiles(self, dataset):
        """
        Get all dummy count filenames for non-split counts
        :param dataset: DROP group name from wildcard
        :return: list of files
        """
        ids = self.sampleAnnotation.getIDsByGroup(dataset, assay="RNA")
        file_stump = self.processedDataDir / "aberrant_splicing" / "datasets" / \
                     "cache" / "nonSplicedCounts" / f"raw-{dataset}" 
        done_files = str(file_stump / "nonSplicedCounts-{sample_id}.h5")
        return expand(done_files, sample_id=ids)


    def getExternalCounts(self, group: str, fileType: str = "k_j_counts"):
        """
        Get externally provided splice count data dir based on the given group.
        If a file type is given the corresponding file within the folder is returned. 
        :param group: DROP group name from wildcard
        :param fileType: name of the file without extension which is to be returned
        :return: list of directories or files
        """
        ids = self.sa.getIDsByGroup(group, assay="SPLICE_COUNT")
        extCountFiles = self.sa.getImportCountFiles(annotation=None, group=group, 
                file_type="SPLICE_COUNTS_DIR", asSet=False)
        if fileType is not None:
            extCountFiles = np.asarray(extCountFiles)[pd.isna(extCountFiles) == False].tolist()
            extCountFiles = [x + "/" + fileType + ".tsv.gz" for x in extCountFiles]
        return extCountFiles
    
