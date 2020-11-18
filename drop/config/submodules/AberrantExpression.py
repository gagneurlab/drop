from snakemake.io import expand
import numpy as np

from drop import utils
from .Submodules import Submodule


class AE(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir)
        self.CONFIG_KEYS = [
            "groups", "fpkmCutoff", "implementation", "padjCutoff", "zScoreCutoff",
            "maxTestedDimensionProportion"
        ]
        self.name = "AberrantExpression"
        self.rnaIDs = self.sa.subsetGroups(self.groups, assay="RNA")
        self.extRnaIDs = self.sa.subsetGroups(self.groups, assay="GENE_COUNTS")

        # check number of IDs per group
        all_ids = {g: self.rnaIDs[g] + self.extRnaIDs[g] for g in self.groups}
        self.checkSubset(all_ids)

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
        extCountFiles = self.sa.getImportCountFiles(annotation, group, file_type="GENE_COUNTS_FILE")
        count_files.extend(extCountFiles)
        return count_files

    def getCountParams(self, rnaID):
        sa_row = self.sa.getRow("RNA_ID", rnaID)
        count_params = sa_row[["STRAND", "COUNT_MODE", "PAIRED_END", "COUNT_OVERLAPS"]]
        count_params_dict = {
            k: bool(v) if isinstance(v, np.bool_) else v
            for k, v in count_params.iloc[0].to_dict().items()
        }
        return count_params_dict
        # count_params.iloc[0].to_dict()