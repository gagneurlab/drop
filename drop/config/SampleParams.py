from pathlib import Path
from snakemake.logging import logger
import numpy as np
import pandas as pd
import os
import shutil



class ParamHelper:
    def __init__(self,include,sampleAnnotationColumns,group,path):
        """
        include: boolean. True-include the sample annotation columns, False- use the compliment columns
        sampleAnnotationColumns: list. List of columns in the sample annotation table to include (or not)
        group: boolean. True- treat the IDs as a drop group, False- treat the IDs separately as individual files
        path: string. path of the directory to write the sample param file to
        """
        self.include= include
        self.sampleAnnotationColumns = sampleAnnotationColumns
        self.group = group
        self.path = path


class SampleParams:
    # dictionary containing the module names, and the predetermined subfolders for the processed data
    MODULE_NAMES = {"AberrantExpression": "aberrant_expression",
                    "AberrantSplicing": "aberrant_splicing",
                    "rnaVariantCalling": "rnaVariantCalling",
                    "MonoallelicExpression": "mae"}


    # helper object containing the relevant information for the module/param pair
    # each Param Helper has [include,SA columns, group,path]
    AE_countParams = ParamHelper(
    # include the following columns (in the list), don't group the samples by DROP/RVC group, and use the suffix _counts
                   True,
                   ["RNA_ID", "RNA_BAM_FILE","COUNT_MODE", "PAIRED_END", "COUNT_OVERLAPS", "STRAND"],
                   False,
                   "counts")

    # include the following columns (in the list), group the samples by DROP/RVC group, and use the suffix _merge
    AE_mergeParams = ParamHelper(
                   True,
                   ["RNA_ID", "RNA_BAM_FILE","GENE_COUNTS_FILE","COUNT_MODE", "PAIRED_END", "COUNT_OVERLAPS", "STRAND","GENE_ANNOTATION"],
                   True,
                   "merge")

    # include the following columns (in the list), group the samples by DROP/RVC group, and use the suffix _results
    AE_resultParams= ParamHelper(
                   True,
                   ["RNA_ID","DNA_ID","HPO_TERMS","GENE_COUNTS_FILE","GENE_ANNOTATION"],
                   True,
                   "results")

    # include the following columns (in the list), don't group the samples by DROP/RVC group, and use the suffix _snvs
    MAE_snvParams = ParamHelper(
                   True,
                   ["RNA_ID","DNA_ID","RNA_BAM_FILE", "DNA_VCF_FILE","GENOME"],
                   False,
                   "snvs")

    # exclude the following columns (in the list), group the samples by DROP/RVC group, and use the suffix _results
    MAE_resultParams = ParamHelper(
                   True,
                   ["RNA_ID","DNA_ID","DROP_GROUP","RNA_BAM_FILE", "DNA_VCF_FILE","GENOME"],
                   True,
                   "results")

    # include the following columns (in the list), don't group the samples by DROP/RVC group, and use the suffix _samples
    RVC_sampleParams = ParamHelper(
                   True,
                   ["RNA_ID","RNA_BAM_FILE","DROP_GROUP","GENOME"],
                   False,
                   "samples")

    # include the following columns (in the list), group the samples by DROP/RVC group, and call it batches
    RVC_batchParams = ParamHelper(
                   True,
                   ["RNA_ID","RNA_BAM_FILE","DROP_GROUP","GENOME"],
                   True,
                   "batches")

    # dictionary containing the key type of parameter and the corresponding param information object
    PARAM_COLS = {"AberrantExpression":
                     { "countParams":  AE_countParams,
                       "mergeParams":  AE_mergeParams,
                       "resultParams": AE_resultParams
                     },

                  "MonoallelicExpression":
                     { "snvParams": MAE_snvParams,
                       "resultParams": MAE_resultParams
                     },

                  "rnaVariantCalling":
                     { "sampleParams": RVC_sampleParams,
                       "batchParams": RVC_batchParams
                     }
                 }


    def __init__(self,AE,AS,MAE,RVC,geneAnnotation,processedDataDir, sampleAnnotation):
        """
        AE: object. AberrantExpression object as created in DropConfig.py
        AS: object. AberrantSplicing object as created in DropConfig.py
        MAE: object. MonoallelicExpression object as created in DropConfig.py
        RVC: object. rnaVariantCalling object as created in DropConfig.py
        annotation: dict. dictionary containing the annotation ID (version) and the path to it
        processedDataDir: string. path to the processedDataDir
        sampleAnnotation: object. SampleAnnotation object as defined by SampleAnnotation.py
        """

        moduleList = [x for x in [AE,AS,MAE,RVC] if x.run]
        self.geneAnnotation = geneAnnotation
        self.processedDataDir = processedDataDir
        self.sampleAnnotation = sampleAnnotation

        # loop through each module and write the sample parameters
        self.writeSAtableSampleParams(moduleList)
        self.writeConfigSampleParams(moduleList)

    def writeConfigSampleParams(self,moduleList):
        """
        moduleList: list. list containing module objects to loop through
        """
        for module in moduleList:
            config_dict = module.dict_
            format_dict = {key:str(config_dict[key]) for key in config_dict}
            current_config = pd.DataFrame.from_dict(format_dict,orient = "index",columns = ["value"])
            moduleCSV = self.processedDataDir / self.MODULE_NAMES[module.name] / "params" / "config" 
            moduleCSV.mkdir(parents = True,exist_ok = True)
            filename = f"{str(moduleCSV)}/{module.name}_config.tsv"

            # if a file by the desired name exists.
            if os.path.isfile(filename):
                old_config = pd.read_csv(filename,sep = "\t",index_col = 0).fillna(value = "").astype(str)
                if current_config.equals(old_config):
                    pass
                else:
                    # if they're different remove the existing file and rename TEMP to the desired file. Updating to the current SA table
                    logger.info(f"{filename} Param Files do not match. Updating to current Sample Annotation\n")
                    current_config.to_csv(filename, sep = "\t", index = True, header = True,na_rep = "")

            else:
                #logger.info("{} Param File did not already exist. Writing it\n".format(filename))
                current_config.to_csv(filename, sep = "\t", index = True, header = True,na_rep = "")



    def writeSAtableSampleParams(self,moduleList):
        """
        moduleList: list. list containing module objects to loop through
        using the dictionary and helpers defined above, loop through each module and param pair and write the sample param file
        """
        for ann in self.geneAnnotation:
            for module in moduleList:
                if module.name == "AberrantExpression":
                    modulePath = self.processedDataDir / self.MODULE_NAMES[module.name] / ann / "params"
                elif module.name == "MonoallelicExpression":
                    modulePath = self.processedDataDir / self.MODULE_NAMES[module.name] / "params"
                elif module.name == "rnaVariantCalling":
                    modulePath = self.processedDataDir / self.MODULE_NAMES[module.name] / "params"
                elif module.name == "AberrantSplicing":
                    continue
                    #modulePath = self.processedDataDir / self.MODULE_NAMES[module.name] / "params"
                else: raise(ValueError,"currently only AberrantExpression, AberrantSplicing, and MonoallelicExpression are supported")

                paramKeys = self.PARAM_COLS[module.name].keys()
                for paramType in paramKeys:
                    self.writeSampleParams(
                             module,
                             modulePath / self.PARAM_COLS[module.name][paramType].path,
                             paramType,
                             self.PARAM_COLS[module.name][paramType].sampleAnnotationColumns,
                             self.PARAM_COLS[module.name][paramType].include,
                             self.PARAM_COLS[module.name][paramType].group,
                         )

    def writeSampleParams(self,module,path,file_suffix,param_cols,include,group_param):
        """
        module: object. Drop module object, used to get drop group attributes
        path: string. path to where to write the param files
        file_suffix: string. what type of sample param are we writing (countParam, resultParam...)
        param_cols: list. list of sample annotation columns to use to determine parameters for that job
        include: boolean. True- include all of the columns in param_cols to build param file. False- use all other columns in SA
        group_param: boolean. True- group by drop_group in sample annotation table. False- treat each sample individually
        """
        # initialize groups and sa table
        module_groups = module.groups
        sa_df = self.sampleAnnotation.annotationTable

        all_RNA_ids = []
        # for each DROP_GROUP used for each module build the list of RNA_IDs that are going to be merged for that run
        # also build a list of all RNA_IDs triggered in this run
        for group in module_groups:
            group_IDs = self.sampleAnnotation.getIDsByGroup(group) + \
                          self.sampleAnnotation.getIDsByGroup(group, assay="GENE_COUNT")
            all_RNA_ids = all_RNA_ids + group_IDs

            if group_param:
                # write the params file for the Merge, and also the info_params file for the Results
                self.updateParamFiles(path,f"{group}_{file_suffix}.csv", \
                                    sa_df,param_cols,group_IDs,include)

            else:
            # for all unique RNA_IDs compiled across the groups build the individual param files used for Counts
                for ID in set(all_RNA_ids):
                    self.updateParamFiles(path,f"{ID}_{file_suffix}.csv", \
                                    sa_df,param_cols,[ID],include)



    def updateParamFiles(self,path,filename,sa_df,param_cols,ID,include):
        """
        path: string. path to where to write the param files (separate from filename to build path if not existing)
        filename: string. name of file to write in path
        param_cols: list. list of sample annotation columns to use to determine parameters for that job
        ID: list. list containing the identifier for the sa_col. either the [sample name] or the [samples in drop group]
        include: boolean. True- include all of the columns in param_cols to build param file. False- use all other columns in SA
        """
        # build the path to the param file
        path.mkdir(parents = True,exist_ok = True)

        param_cols = [col for col in param_cols if col in sa_df.columns] # remove params that are not columns in SA table

        # take the complement of columns if indicated by !include
        if not include:
            param_cols = [col for col in sa_df.columns if col not in param_cols]

        # designate the TEMP and final param file names
        true_filename = f"{path}/{filename}"


        # if a file by the desired name exists.
        if os.path.isfile(true_filename):

            # replace any strings of nan with "NA"
            current_SA = sa_df.loc[sa_df["RNA_ID"].isin(ID),param_cols].reset_index(drop = True)
            current_SA = current_SA.replace("nan","NA").fillna(value = "NA").astype(str)
            old_SA = pd.read_csv(true_filename).reset_index(drop = True).fillna(value = "NA").astype(str)

            if current_SA.equals(old_SA):
                pass
            else:
                # if they're different remove the existing file and rename TEMP to the desired file. Updating to the current SA table
                logger.info(f"{filename} Param Files do not match. Updating to current Sample Annotation\n")
                current_SA.to_csv(true_filename, index = False,header = True,na_rep = "NA")
        # if the param file doesn't exist, just write to the desired file
        else:
            #logger.info("{} Param File did not already exist. Writing it\n".format(filename))
            sa_df.loc[sa_df["RNA_ID"].isin(ID),param_cols].to_csv(true_filename, index = False,header = True,na_rep = "NA")
