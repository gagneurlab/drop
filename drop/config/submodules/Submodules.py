from pathlib import Path
from drop import utils
from snakemake.logging import logger


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

    def checkSubset(self, groupSubsets, warn=30, error=10):
        """
        Give warning or error if subsetting results in too few sample IDs per group.
        :param groupSubsets:
        :param warn: number of samples threshold at which to warn about too few samples
        :param error: number of samples threshold at which to give error
        """
        for group in self.groups:
            if len(groupSubsets[group]) < error:
                message = f'Too few IDs in DROP_GROUP {group}'
                message += f', please ensure that it has at least {error} IDs'
                message += f', groups: {groupSubsets[group]}'
                raise ValueError(message)
            elif len(groupSubsets[group]) < warn:
                logger.info(f'WARNING: Less than {warn} IDs in DROP_GROUP {group}')
