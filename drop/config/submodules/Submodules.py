from pathlib import Path
from drop import utils


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

