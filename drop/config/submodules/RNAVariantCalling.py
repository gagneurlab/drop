from drop import utils
from .Submodules import Submodule

class RVC(Submodule):

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir,workDir,genome):
        super().__init__(config, sampleAnnotation, processedDataDir, processedResultsDir,workDir)
        self.CONFIG_KEYS = [
            "run","groups","KGsnps", "millsIndels", "dbSNP","repeat_mask",
            "hcArgs","maxVarFreqCohort","createSingleVCF","minAlt"
        ]
        self.name = "rnaVariantCalling"
        # if self.run is false return without doing any config/sa checks for completeness
        if not self.run:
            return
        self.rnaIDs = self.sampleAnnotation.subsetGroups(self.groups)
        self.batchIDs = self.setBatchDict()
        self.batch_genome = {}
        # genomeFiles{config_name -> path} from config and sampleGenomes {sampleID -> config_name} from SA
        self.genomeFiles = genome.reference
        self.sampleGenomes = self.setGenomeDict(self.genomeFiles)

        #update sampleGenomes dict with key for the batch itself
        self.sampleGenomes.update(self.check_batch_genome())

    def check_batch_genome(self):
        """
        raise an error if a batch as definied by the sample annotation column contains more than one genome
        as defined by the sample annotation column GENOME
        """
        ref_genomes = dict()
        for batch in self.sampleAnnotation.rnaIDs:
            if batch in self.groups:
                ref_genomes[batch] = set()
                for sample in self.sampleAnnotation.rnaIDs[batch]:
                    ref_genomes[batch].add(self.sampleGenomes[sample])
                if len(ref_genomes[batch]) != 1:
                    raise ValueError(f"The reference genome files {ref_genomes} within a batch must be the same")
                else:
                    ref_genomes[batch] = list(ref_genomes[batch])[0]

        return ref_genomes

    def setBatchDict(self):
        """
        build a dictionary that connects the RNAID to the batch value. sample1 -> batch0 as defined by the sample annotation table
        each sample can only belong to a single batch. sample1 -> batch0 sample1 -> batch1
        """
        if not self.run:
            return {}
        if (len(self.rnaIDs.items()) < 1): raise ValueError("No RNA IDs found in the group, can not create dictionary")
        dict_ = super().setDefaultKeys(None)
        setKey = utils.setKey
        for key,values in self.rnaIDs.items():
            for v in values:
                if(v) in dict_.keys(): raise ValueError("RNA IDs must be unique to the RNA Variant Calling group, can not have an RNA_ID point to multiple DROP_GROUP groups")
                setKey(dict_, None, v ,key)
        return dict_

    def setDefaultKeys(self, dict_):
        super().setDefaultKeys(dict_)
        setKey = utils.setKey
        setKey(dict_, None, "run", False)
        setKey(dict_, None, "groups", self.sampleAnnotation.getGroups())
        setKey(dict_, None, "highQualityVCFs", [])
        setKey(dict_, None, "dbSNP","") 
        setKey(dict_, None, "repeat_mask", "")
        setKey(dict_, None, "hcArgs", "")
        setKey(dict_, None, "addAF", False)
        setKey(dict_, None, "maxAF", 0.001 )
        setKey(dict_, None, "maxVarFreqCohort", 0.05 )
        setKey(dict_, None, "minAlt", 3)
        setKey(dict_, None, "createSingleVCF", False)
        setKey(dict_, None, "yieldSize", 100000)

        if dict_["run"]:
            dict_ = utils.checkKeys(dict_, keys=["repeat_mask","highQualityVCFs","dbSNP"], check_files=True)
        return dict_
