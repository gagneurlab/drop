from drop import utils
from .Submodules import Submodule

class RVC(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir,workDir,genome):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir,workDir)
        self.CONFIG_KEYS = [
            "groups","KGsnps", "millsIndels", "dbSNP","repeat_mask","hcArgs","minAlt"
        ]
        self.name = "rnaVariantCalling"
        self.rnaIDs = self.sampleAnnotation.subsetGroups(self.groups, assay="RVC")
        self.batchIDs = self.setBatchDict()
        self.batch_genome = {}
        # genomeFiles{config_name -> path} from config and sampleGenomes {sampleID -> config_name} from SA
        self.genomeFiles = genome.reference
        self.sampleGenomes = self.setGenomeDict(self.genomeFiles,group_key="RNA_VARIANT_GROUP")
 
        #update sampleGenomes dict with key for the batch itself
        self.sampleGenomes.update(self.check_batch_genome())


    # check each batch for one reference file per batch
    def check_batch_genome(self):
        ref_genomes = dict()
        for batch in self.sampleAnnotation.rnaIDs_RVC:
            if batch in self.groups:
                ref_genomes[batch] = set()
                for sample in self.sampleAnnotation.rnaIDs_RVC[batch]:
                    ref_genomes[batch].add(self.sampleGenomes[sample])
                if len(ref_genomes[batch]) != 1:
                    raise ValueError(f"The reference genome files {ref_genomes} within a batch must be the same")
                else:
                    ref_genomes[batch] = list(ref_genomes[batch])[0]

        return ref_genomes

    def setBatchDict(self):
        if not self.run:
            return {}
        if (len(self.rnaIDs.items()) < 1): raise ValueError("No RNA IDs found in the group, can not create dictionary")
        dict_ = super().setDefaultKeys(None)
        setKey = utils.setKey
        for key,values in self.rnaIDs.items():
            for v in values:
                if(v) in dict_.keys(): raise ValueError("RNA IDs must be unique to the RNA Variant Calling group, can not have an RNA_ID point to multiple RNA_VARIANT_GROUPs")
                setKey(dict_, None, v ,key)
        return dict_

    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        dict_ = utils.checkKeys(dict_, keys=["repeat_mask"], check_files=True)
        groups = setKey(dict_, None, "groups", self.sampleAnnotation.getGroups(assay="RVC"))
        setKey(dict_, None, "knownVCFs", [])
        setKey(dict_, None, "repeat_mask", "")
        setKey(dict_, None, "hcArgs", "")
        setKey(dict_, None, "minAlt", 3)
        return dict_

