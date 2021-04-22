from drop import utils
from .Submodules import Submodule
from snakemake import logger


class MAE(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir,genomeFiles):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir)
        self.CONFIG_KEYS = [
            "groups", "genome", "qcVcf", "qcGroups", "gatkIgnoreHeaderCheck", "padjCutoff",
            "allelicRatioCutoff", "maxAF", "gnomAD"
        ]
        self.name = "MonoallelicExpression"
        self.qcGroups = self.dict_["qcGroups"]
        self.qcVcfFile = self.dict_["qcVcf"]
        self.maeIDs = self.createMaeIDS(id_sep='--')
        
        # genomeFiles{config_name -> path} from config and sampleGenomes {sampleID -> config_name} from SA
        self.genomeFiles = self.setGenomeFile(genomeFiles)
        self.sampleGenomes = self.setGenomeDict(self.genomeFiles)
        if len(set(self.genomeFiles.keys()) - set(self.sampleGenomes.values())) > 0: 
            logger.error("The genome keys defined in the config do not match exactly the values in the GENOME column of the sample annotation. Please fix them.")
        elif "GENOME" not in self.sa.sa.columns.values: #genome column not defined 
            pass
        elif not all(self.sa.sa["GENOME"] == "nan") and len(self.genomeFiles) == 1:
            logger.error("WARNING: The genome is defined globally in the config, however non-empty values are in the sample annotation table. Using the globally defined path, please consider fixing this.")

    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        dict_ = utils.checkKeys(dict_, keys=["qcVcf"], check_files=True)
        groups = setKey(dict_, None, "groups", self.sa.getGroups(assay="DNA"))
        setKey(dict_, None, "qcGroups", groups)
        setKey(dict_, None, "gatkIgnoreHeaderCheck", True)
        setKey(dict_, None, "padjCutoff", .05)
        setKey(dict_, None, "allelicRatioCutoff", 0.8)
        setKey(dict_, None, "maxAF", .001)
        setKey(dict_, None, "addAF", False)
        setKey(dict_, None, "maxVarFreqCohort", 0.04)
        setKey(dict_, None, "gnomAD", False)
        return dict_

    def createMaeIDS(self, id_sep='--'):
        """
        Create MAE IDs from sample annotation
        :param id_sep: separator
        :return: {drop group name : list of MAE IDs per group}
        """
        grouped_rna_ids = self.sa.subsetGroups(self.groups, assay="RNA")
        self.checkSubset(grouped_rna_ids, warn=1, error=1)
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
        Useful for collecting all MAE IDs ungrouped.
        :return: list of all MAE IDs from the groups specified in the config
        """
        all_ids = []
        groups = [self.groups] if isinstance(self.groups, str) else self.groups
        for group in groups:
            all_ids.extend(self.getMaeByGroup(group))
        return all_ids

    def getVcf(self, id):
        if id == 'QC':
            return self.qcVcfFile
        return self.sa.getFilePath(id, 'DNA_VCF_FILE')

     
    # set the genomeFiles to be a dictionary. if globally (or MAE) defined as a string. force into dictionary {genome:genome}
    def setGenomeFile(self,genomeFiles):
        if isinstance(genomeFiles,str):
            return {genomeFiles:genomeFiles}
        else:
            return genomeFiles

    # map out the samples in the group to the corresponding genome defined in SA
    def setGenomeDict(self,genomeFiles):
        genomeDict = {}
        if len(genomeFiles) == 1: #globally defined in the config
            globalGenome = list(genomeFiles.values())[0]

            # subset SA by the drop group (not exact match) and skip the filtering by SA-GENOME column
            genomeDict = self.sa.getGenomes(globalGenome,self.groups, file_type="RNA_ID",
                                            column="DROP_GROUP", group_key="DROP_GROUP",exact_match = False,skip = True)
        else:
            # subset SA by the drop group (not exact match) and filter by SA-GENOME column. Must exactly match config key
            for gf in genomeFiles.keys():
                genomeDict.update(self.sa.getGenomes(gf,self.groups, file_type="RNA_ID",
                                            column="GENOME", group_key="DROP_GROUP",exact_match = False,skip = False))

        return genomeDict  

    # look up for a sampleID genomeFiles{ncbi -> path} and sampleGenomes {sampleID -> ncbi}
    def getGenomePath(self,sampleID):
        try:
            return self.genomeFiles[self.sampleGenomes[sampleID]]
        except KeyError:
            raise KeyError(f"The Config file has defined specific key,value for genome path \
but the SA table does not match for sample {sampleID}")
