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
        """
        Rename the local Scripts directory based on whether or not a module has been designated to run in the config
        wBuild ignores folders with the '_' prefix, so add those to all Script files that won't be run. For modules
        that will be turned on now that were previously off, remove the '_' prefix
        """
        work_dir_prefix = utils.returnPath(self.workDir / "Scripts", False)
        work_dir_off = work_dir_prefix / ("_" + self.name)
        work_dir_on = work_dir_prefix / self.name

        # if the module will be run, rename to the {work_dir_on} naming
        if self.run:
            work_dir = Path(work_dir_on / "pipeline")
            if (os.path.isdir(work_dir_off)):
                try:
                    os.rename(work_dir_off, work_dir_on)
                except:
                    raise OSError(f"Could not rename {work_dir_off} to {work_dir_on}\n \
                        try running 'drop update' to reset file structure\n")

        # if the module will be not run, rename to the {work_dir_off} naming
        else:
            logger.info(f"{self.name} has been turned off in the config file")
            work_dir = Path(work_dir_off / "pipeline")
            if (os.path.isdir(work_dir_on)):
                try:
                    os.rename(work_dir_on, work_dir_off)
                except:
                    raise OSError(f"Could not rename {work_dir_on} to {work_dir_off}\n \
                        try running 'drop update' to reset file structure\n")

        return work_dir

    def getModuleIndexFiles(self, index_name, work_dir):
        """
        Replace the rule in the Snakefile to be run by each module. Define the IndexFile for each module
        :param index_name: name of the html index file
        :param work_dir: the working dir as defined in the config
        """
        if (self.run):
            index_input, index_output, graph_file, _ = createIndexRule(
                scriptsPath=str(work_dir),
                index_name=index_name
            )
            return index_input, graph_file, index_output
        else:
            # dummy return of hidden pipeline directory which already exists
            return self.getWorkdir(hide_dir=True), self.getWorkdir(hide_dir=True), self.getWorkdir(hide_dir=True)


    def setGenomeDict(self, genomeFiles,group_key="DROP_GROUP"):
        """
        create a dictionary that connects the module groups to the approrpiate genome based on the sample annotation table
        :param genomeFiles: the genome files as defined in the config 'genome'. Either a dictionary or a string
        :param group_key: by default "DROP_GROUP" this defines the groups within the group_key column
        """
        genomeDict = {}
        if len(genomeFiles) == 1:  # globally defined in the config
            globalGenome = list(genomeFiles.values())[0]

            # subset SA by the drop group and skip the filtering by SA-GENOME column
            # because only 1 genome is defined don't filter by the GENOME column (it may not exist)
            genomeDict = self.sampleAnnotation.getGenomes(
                globalGenome,
                self.groups,
                file_type="RNA_ID",
                column=group_key, group_key=group_key,
                skip=True
            )
        else:
            # subset SA by the drop group and filter by SA-GENOME column. Must exactly match config key
            # because more than 1 genome is defined filter by the GENOME column 
            for gf in genomeFiles.keys():
                genomeDict.update(
                    self.sampleAnnotation.getGenomes(
                        gf,
                        self.groups,
                        file_type="RNA_ID",
                        column="GENOME", group_key=group_key,
                        skip=False
                    )
                )

        return genomeDict

    def getGenomePath(self, sampleID):
        """
        create a dictionary that connects the sample to the approrpiate genome as named in the config
        :param sampleID: the sampleID
        """

        # if only 1 genome is defined, all samples must use that genome. Return the value
        if len(self.genomeFiles) ==1:
            return(list(self.genomeFiles.values())[0])
        # if more than 1 genome is defined use the sample annotation table to match sample an genome
        else:
            try:
                return self.genomeFiles[self.sampleGenomes[sampleID]]
            except KeyError:
                raise KeyError(
                    f"The Config file has defined specific key,value for genome path "
                    f"but the SA table does not match for sample {sampleID}"
                )
