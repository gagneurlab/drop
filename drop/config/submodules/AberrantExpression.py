from snakemake.io import expand
import numpy as np

from drop import utils
from .Submodules import Submodule


class AE(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir, workDir):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir, workDir)
        self.CONFIG_KEYS = [
            "groups", "fpkmCutoff", "implementation", "padjCutoff", "zScoreCutoff",
            "maxTestedDimensionProportion"
        ]
        self.name = "AberrantExpression"
        # if self.run is false return without doing any config/sa checks for completeness
        if not self.run:
            return
        self.rnaIDs = self.sampleAnnotation.subsetGroups(self.groups, assay="RNA")
        self.extRnaIDs = self.sampleAnnotation.subsetGroups(self.groups, assay="GENE_COUNTS")
        for g in self.groups:
            if len(set(self.rnaIDs[g]) & set(self.extRnaIDs[g])) > 0:
                raise ValueError(f"{set(self.rnaIDs[g]) & set(self.extRnaIDs[g])} has both BAM and external count file \
                please fix to only have either external count or BAM processing\n")

        # check number of IDs per group
        all_ids = self.sampleAnnotation.subsetGroups(self.groups, assay=["RNA", "GENE_COUNTS"])
        self.checkSubset(all_ids)

    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        setKey(dict_, None, "run", False)
        setKey(dict_, None, "groups", self.sampleAnnotation.getGroups(assay="RNA"))
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


        bam_IDs = self.sampleAnnotation.getIDsByGroup(group, assay="RNA")
        file_stump = self.processedDataDir / "aberrant_expression" / annotation / "counts" / "{sampleID}.Rds"
        count_files = expand(str(file_stump), sampleID=bam_IDs)
        # if sample annotation table does not contain GENE_COUNTS_FILE column. return no external counts
        if("GENE_COUNTS_FILE" not in self.sampleAnnotation.SAMPLE_ANNOTATION_COLUMNS):
            extCountFiles = []
        else:
            extCountFiles = self.sampleAnnotation.getImportCountFiles(annotation, group, file_type="GENE_COUNTS_FILE")
            count_files.extend(extCountFiles)
        return count_files

    def getCountParams(self, rnaID):
        sa_row = self.sampleAnnotation.getRow("RNA_ID", rnaID)
        count_params = sa_row[["STRAND", "COUNT_MODE", "PAIRED_END", "COUNT_OVERLAPS"]]
        count_params_dict = {
            k: bool(v) if isinstance(v, np.bool_) else v
            for k, v in count_params.iloc[0].to_dict().items()
        }
        return count_params_dict
