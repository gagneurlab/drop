from drop import utils
from .Submodules import Submodule


class MAE(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir)
        self.CONFIG_KEYS = [
            "groups", "genome", "qcVcf", "qcGroups", "gatkIgnoreHeaderCheck", "padjCutoff",
            "allelicRatioCutoff", "maxAF", "gnomAD"
        ]
        self.name = "MonoallelicExpression"
        self.qcGroups = self.dict_["qcGroups"]
        self.qcVcfFile = self.dict_["qcVcf"]
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

    def getQcVcf(self):
        return self.qcVcfFile
