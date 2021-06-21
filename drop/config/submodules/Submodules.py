from drop import utils
from snakemake.logging import logger
import os
from pathlib import Path
from wbuild.createIndex import createIndexRule


class Submodule:

    def __init__(self, config, sampleAnnotation, processedDataDir, processedResultsDir, workDir):
        self.CONFIG_KEYS = []
        self.name = "Submodule"
        self.sampleAnnotation = sampleAnnotation
        self.processedDataDir = processedDataDir
        self.processedResultsDir = processedResultsDir
        self.workDir = workDir
        self.dict_ = self.setDefaultKeys(config)
        self.groups = self.dict_["groups"]
        self.run = self.dict_.get('run', True)

    def setDefaultKeys(self, dict_):
        dict_ = {} if dict_ is None else dict_
        return dict_

    def get(self, key):
        if key not in self.CONFIG_KEYS:
            raise KeyError(f"{key} not defined for {self.name} config")
        return self.dict_[key]

    def getWorkdir(self, hide_dir=False, str_=True):
        if hide_dir:
            return utils.returnPath(self.workDir / "Scripts" / ("_" + self.name) / "pipeline", str_)
        else:
            return utils.returnPath(self.workDir / "Scripts" / self.name / "pipeline", str_)

    def getSnakefile(self):
        if os.path.isdir(self.getWorkdir()):
            return self.getWorkdir() + "/Snakefile"
        else:
            return self.getWorkdir(hide_dir=True) + "/Snakefile"

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

    def renameLocalDir(self):
        work_dir_prefix = utils.returnPath(self.workDir / "Scripts", False)
        work_dir_off = work_dir_prefix / ("_" + self.name)
        work_dir_on = work_dir_prefix / self.name

        if self.run:
            work_dir = Path(work_dir_on / "pipeline")
            if (os.path.isdir(work_dir_off)):
                try:
                    os.rename(work_dir_off, work_dir_on)
                except:
                    raise OSError(f"Could not rename {work_dir_off} to {work_dir_on}\n \
                        try running 'drop update' to reset file structure\n")
        else:
            work_dir = Path(work_dir_off / "pipeline")
            if (os.path.isdir(work_dir_on)):
                try:
                    os.rename(work_dir_on, work_dir_off)
                except:
                    raise OSError(f"Could not rename {work_dir_on} to {work_dir_off}\n \
                        try running 'drop update' to reset file structure\n")

        return work_dir

    def getModuleIndexFiles(self, index_name, work_dir):
        run = self.run
        if (run):
            index_input, index_output, graph_file, _ = createIndexRule(
                scriptsPath=str(work_dir),
                index_name=index_name
            )
            return index_input, graph_file, index_output
        else:
            # dummy return of hidden pipeline directory which already exists
            return self.getWorkdir(hide_dir=True), self.getWorkdir(hide_dir=True), self.getWorkdir(hide_dir=True)
