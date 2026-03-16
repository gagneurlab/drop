from snakemake.io import expand
from snakemake.logging import logger
import numpy as np

from drop import utils
from .Submodules import Submodule

class AE(Submodule):

    def __init__(self, config, sampleAnnotation: "SampleAnnotation", processedDataDir, processedResultsDir, workDir):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir, workDir)
        self.CONFIG_KEYS = [
            "groups", "fpkmCutoff", "implementation", "padjCutoff", "zScoreCutoff",
            "maxTestedDimensionProportion", "genesToTest", "useOHTtoObtainQ"
        ]
        self.name = "AberrantExpression"
        # if self.run is false return without doing any config/sa checks for completeness
        if not self.run:
            return
        self.rnaIDs = self.sampleAnnotation.subsetGroups(self.groups, assay="RNA")
        self.extRnaIDs = self.sampleAnnotation.subsetGroups(self.groups, assay="GENE_COUNTS")
        # Note: mixing BAM and external counts is now allowed to support using external expression counts
        # while still processing BAM files for splicing analysis

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
        setKey(dict_, None, "genesToTest", None)
        setKey(dict_, None, "maxTestedDimensionProportion", 3)
        setKey(dict_, None, "yieldSize", 2000000)
        setKey(dict_, None, "useOHTtoObtainQ", True)
        return dict_

    def getCountFiles(self, annotation, group):
        """
        Get all count files from DROP (counted from BAM file) and external count matrices
        When both BAM and external counts are available, prioritize external counts for expression analysis
        while still allowing BAM files to be used for splicing analysis
        :param annotation: annotation name from wildcard
        :param group: DROP group name from wildcard
        :return: list of files
        """

        bam_IDs = self.sampleAnnotation.getIDsByGroup(group, assay="RNA")
        ext_IDs = self.sampleAnnotation.getIDsByGroup(group, assay="GENE_COUNTS")
        
        # Get BAM-based count files only for samples that don't have external counts
        bam_only_IDs = [id_ for id_ in bam_IDs if id_ not in ext_IDs]
        
        # Log information about mixed usage when both BAM and external counts are available
        if ext_IDs and bam_only_IDs:
            logger.info(f"Group '{group}': Using external expression counts for {len(ext_IDs)} samples "
                       f"({', '.join(ext_IDs)}) and BAM-derived counts for {len(bam_only_IDs)} samples "
                       f"({', '.join(bam_only_IDs)})")
        elif ext_IDs:
            logger.info(f"Group '{group}': Using external expression counts for all {len(ext_IDs)} samples")
        
        file_stump = self.processedDataDir / "aberrant_expression" / annotation / "counts" / "{sampleID}.Rds"
        count_files = expand(str(file_stump), sampleID=bam_only_IDs)
        
        # Add external count files if available
        if("GENE_COUNTS_FILE" in self.sampleAnnotation.SAMPLE_ANNOTATION_COLUMNS):
            extCountFiles = self.sampleAnnotation.getImportCountFiles(annotation, group, file_type="GENE_COUNTS_FILE")
            count_files.extend(extCountFiles)
            
        return count_files

    def getCountParams(self, rnaID) -> dict[str, str | bool]:
        sa_row = self.sampleAnnotation.getRow("RNA_ID", rnaID)
        return {
            "STRAND": sa_row.iloc[0].STRAND,
            "COUNT_MODE": sa_row.iloc[0].COUNT_MODE,
            "PAIRED_END": self._convert_to_bool(sa_row.iloc[0].PAIRED_END),
            "COUNT_OVERLAPS": self._convert_to_bool(sa_row.iloc[0].COUNT_OVERLAPS),
        }

    @staticmethod
    def _convert_to_bool(v: object) -> bool:
        """Convert value to bool if possible, otherwise return the original object."""
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        if isinstance(v, (int, float)):
            if v == 1:
                return True
            if v == 0:
                return False
        if isinstance(v, str) and v.lower() in ('true', 'false', 'yes', 'no', 't', 'f', 'y', 'n', '1', '0'):
            return v.lower() in ('true', 'yes', 't', 'y', '1')
        raise ValueError(f"Cannot convert value '{v}' of type {type(v)} to bool." +
                         " Please provide a boolean value (True/False), 1/0 or yes/no.")

# import in the end to overcome circular import issues
from drop.config.SampleAnnotation import SampleAnnotation