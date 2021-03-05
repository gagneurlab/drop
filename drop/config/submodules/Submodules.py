from pathlib import Path
from drop import utils
from snakemake.logging import logger
import numpy as np                                                                                                      
import pandas as pd
import os                                                                                                               
import shutil


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

    def update_param_files(self,param_path,filename,sa_df,param_cols,ID,sa_col = "RNA_ID",include = True):
        # build the path to the param file
        param_path.mkdir(parents = True,exist_ok = True)

        param_cols = [col for col in param_cols if col in sa_df.columns] # remove params that are not columns in SA table

        # take the complement of columns if indicated by !include
        if not include:
            param_cols = [col for col in sa_df.columns if col not in param_cols]
    
        # designate the TEMP and final param file names
        true_filename = "{param_path}/{filename}".format(param_path = param_path, filename = filename)

    
        # if a file by the desired name exists. 
        if os.path.isfile(true_filename):

            # replace any strings of nan with "NA"
            current_SA = sa_df.loc[sa_df[sa_col].isin(ID),param_cols].reset_index(drop = True)
            current_SA = current_SA.replace("nan","NA").fillna(value = "NA")
            old_SA = pd.read_csv(true_filename).reset_index(drop = True).fillna(value = "NA")

            if current_SA.equals(old_SA):
                pass
            else:
                # if they're different remove the existing file and rename TEMP to the desired file. Updating to the current SA table
                logger.info("{} Param Files do not match. Updating to current Sample Annotation\n".format(filename))
                current_SA.to_csv(true_filename, index = False,header = True,na_rep = "NA")
                
        # if the param file doesn't exist, just write to the desired file
        else:
            logger.info("{} Param File did not already exist. Writing it\n".format(filename))
            sa_df.loc[sa_df[sa_col].isin(ID),param_cols].to_csv(true_filename, index = False,header = True,na_rep = "NA")

    def writeSampleParams(self,path,params,include,file_suffix,group_param = False):
        # initialize groups and sa table                                                                                
        module_groups = self.groups                                                                                     
        sa_df = self.sa.sa                                                                                              
                                                                                                                        
        all_RNA_ids = []                                                                                                
        # for each DROP_GROUP used for aberrantExpression build the list of RNA_IDs that are going to be merged for that run
        # also build a list of all RNA_IDs triggered in this run                                                        
        for group in module_groups:                                                                                     
            group_IDs = self.sa.getIDsByGroup(group, assay="RNA") + \
                          self.sa.getIDsByGroup(group, assay="GENE_COUNT")                                              
            all_RNA_ids = all_RNA_ids + group_IDs                                                                       
                                                                                                                        
            if group_param:                                                                                             
                # write the params file for the Merge, and also the info_params file for the Results                    
                self.update_param_files(path,f"{group}_{file_suffix}.csv", \
                                    sa_df,params,group_IDs,sa_col = "RNA_ID",include = include)                         
                                                                                                                        
            else:                                                                                                       
            # for all unique RNA_IDs compiled across the groups build the individual param files used for Counts        
                for ID in set(all_RNA_ids):                                                                             
                    self.update_param_files(path,f"{ID}_{file_suffix}.csv", \
                                    sa_df,params,[ID],sa_col = "RNA_ID",include = include)                        
                                                                                                   
