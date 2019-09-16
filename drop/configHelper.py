import os
import pandas as pd
import wbuild.utils as wbu
from snakemake.logging import logger
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')

class ConfigHelper:
    
    def __init__(self, config, html_root=None):

        if not config:
            wconf = wbu.Config()
            config = wconf.conf_dict
        if html_root:
            config["htmlOutputPath"] = f"{html_root}/{config['htmlOutputPath']}"
        
        # set default parameters for missing keys
        if "indexWithFolderName" not in config:
            config["indexWithFolderName"] = True
        if "fileRegex" not in config:
            config["fileRegex"] = ".*\.R"
        if "use_gene_names" not in config:
            config["use_gene_names"] = True
        
        self.config = config
        
        # SAMPLE_FILE_MAPPING has to have the following structure:
        #   [ID | FILE | ASSAY ] , ASSAY can be for example 'RNA_Seq'
        df_mapping = pd.read_csv(config["SAMPLE_FILE_MAPPING"], sep='\t')
        if not set(df_mapping.columns.values)=={"ID", "FILE", "ASSAY"}:
            raise ValueError(f"File columns {df_mapping.columns.values} of sample-file mapping do not correspond to " +
            f"required format with columns [ID | FILE | ASSAY]. File: {config['SAMPLE_FILE_MAPPING']}")
        
        # cleaning SAMPLE_FILE_MAPPING
        df_mapping.dropna(inplace=True)
        df_mapping["existent"] = [os.path.exists(x) for x in df_mapping["FILE"]]
        df_mapping = df_mapping[df_mapping["existent"]]
        df_mapping.drop_duplicates(inplace=True)
        
        self.sample_file_mapping = df_mapping
        
        # SAMPLE_ANNOTATION must have assay names as specified in sample-file mappping for ID columns
        sa_file = config["SAMPLE_ANNOTATION"]
        logger.debug(f"Loading annotation file: '{sa_file}'")
        self.sample_annotation = pd.read_csv(sa_file, sep='\t')
        
        # Group IDs
        # remove unwanted characters
        self.sample_annotation["DROP_GROUP"] = self.sample_annotation["DROP_GROUP"].str.replace("(", "").str.replace(")", "")
        self.all_rna_ids = self.createGroupIds(group_key="DROP_GROUP", assay_key="RNA_ASSAY", sep=',')
        
        ## outrider
        if "outrider_groups" not in config:
            self.config["outrider_groups"] = list(self.all_rna_ids.keys())
        if "min_outrider_ids" not in config or config["min_outrider_ids"] is None:
            self.config["min_outrider_ids"] = 40
        self.outrider_all = self.subsetGroups(self.all_rna_ids, self.config["outrider_groups"])
        self.outrider_filtered = {name:ids for name, ids in self.outrider_all.items() if len(ids) > self.config["min_outrider_ids"]}
        config["outrider_all"], config["outrider_filtered"] = self.outrider_all, self.outrider_filtered
        
        ## fraser
        if "fraser_groups" not in config:
            self.config["fraser_groups"] = list(self.all_rna_ids.keys())
        if "min_fraser_ids" not in config or config["min_fraser_ids"] is None:
            self.config["min_fraser_ids"] = 40
        self.fraser_all = self.subsetGroups(self.all_rna_ids, self.config["fraser_groups"])
        self.fraser_filtered = {name:ids for name, ids in self.fraser_all.items() if len(ids) > self.config["min_fraser_ids"]}
        config["fraser_all"], config["fraser_filtered"] = self.fraser_all, self.fraser_filtered
        
        ## mae
        if "mae_groups" not in config:
            self.config["mae_groups"] = list(self.all_rna_ids.keys())
        mae_rna_by_group = self.subsetGroups(self.all_rna_ids, self.config["mae_groups"])
        self.mae_ids = self.createMaeIDS(mae_rna_by_group, id_sep='--')
        config["mae_ids"] = self.mae_ids
        
    
    """ 
    Get directory path for processed data
    """
    def getProcDataDir(self):
        return self.config["ROOT"] + "/" + self.config["DATASET_NAME"] + "/processed_data"
    
    """ 
    Get directory path for processed results
    """
    def getProcResultsDir(self):
        return self.config["ROOT"] + "/" + self.config["DATASET_NAME"] + "/processed_results"
    
    """
    Get sample ID by experiment
    """
    def getSampleIDs(self, experiment):
        return list(self.sample_file_mapping[self.sample_file_mapping["ASSAY"] == experiment]["ID"]) 
    
    
    def checkFileExists(self, sampleID, assay, verbose=True):
        # note: we already checked for non-existing files in the init, so we only need to check whether the ID is in the sample_file_mapping here
        x = self.sample_file_mapping.query("(ASSAY == @assay) & (ID == @sampleID)")["FILE"]
        exists = (len(x) != 0)
        if (not exists) and verbose:
            print(f"FILE NOT FOUND FOR sampleID: {sampleID} and assay {assay}")
        return exists
        
    """
    Returns vcf and rna files for MAE pipeline
    """
    def getMaeByGroup(self, group):
        if not isinstance(group, str):
            group = list(group)[0]
        return self.mae_ids[group]
    
    def createMaeIDS(self, rna_id_by_group, id_sep='--'):
        
        all_mae_files = self.allMaeFiles()
        
        # subset by group
        mae_ids = {}
        for gr, rna_ids in rna_id_by_group.items():
            mae_subset = all_mae_files [all_mae_files ["RNA_ASSAY"].isin(rna_ids)]
            vcf_rna_pairs = zip(mae_subset["DNA_ID"], mae_subset["RNA_ASSAY"])
            mae_ids[gr] = list(map(id_sep.join, vcf_rna_pairs))
        
        return mae_ids
        
    def allMaeFiles(self):
        # assay ID column in sample_annotation, entry for ASSAY in sample_file_mapping
        
        ####### PROBLEM IF WGS_ASSAY NOT IN SAMPLE_ANNOTATION: possible solution: create empty cols
        if "RNA_ASSAY" not in self.sample_annotation:
            self.sample_annotation["RNA_ASSAY"] = ""
        if "WES_ASSAY" not in self.sample_annotation:
            self.sample_annotation["WES_ASSAY"] = ""
        if "WGS_ASSAY" not in self.sample_annotation:
            self.sample_annotation["WGS_ASSAY"] = ""
        
        mae_files = self.sample_annotation[["RNA_ASSAY", "WES_ASSAY", "WGS_ASSAY"]]
        mae_files = pd.melt(mae_files, id_vars=["RNA_ASSAY"], value_vars=["WES_ASSAY", "WGS_ASSAY"], var_name="DNA_assay", value_name="DNA_ID")
        mae_files.dropna(inplace=True)
        
        # remove IDs of non-existing files
        mae_files['vcf_exists'] = [self.checkFileExists(row["DNA_ID"], row["DNA_assay"], verbose=False) for index, row in mae_files.iterrows()]
        mae_files['rna_exists'] = [self.checkFileExists(x, "RNA_ASSAY") for x in mae_files["RNA_ASSAY"]]
        mae_files = mae_files.query("vcf_exists & rna_exists")
        
        return mae_files
    
    def getRNAByGroup(self, group):
        if not isinstance(group, str):
            group = list(group)[0]
        return self.all_rna_ids[group]
    
    def getFilePath(self, sampleId, assay):
        """
        Function for getting the file path given the sampleId and assay
        @param sampleId: ID of sample
        @param assay: either "RNA_ASSAY", "dna_assay", as specified in the config
        """
        if isinstance(assay, str):
            path = self.sample_file_mapping.query("ASSAY == @assay")
        else:
            path = self.sample_file_mapping[self.sample_file_mapping["ASSAY"].isin(assay)]
        path = path.query("ID == @sampleId")["FILE"]
        return path.iloc[0]
    
    def getFilePaths(self, assay, group=None, ids_by_group=None):
        """
        assay: name of ass column
        group: name of dataset/ group
        ids_by_group: dictionary of IDs by group names
        """
        
        # subset by group if group is specified
        if group is None or ids_by_group is None:
            sampleIDs = self.sample_file_mapping.query("ASSAY == @assay")["ID"]
        else:
            sampleIDs = ids_by_group[group]
        
        files = [] # collect file names
        for sampleID in sampleIDs:
            files.append(self.getFilePath(sampleID, assay))
        return files
    
    """
    Create a full and filtered list of RNA assay IDs subsetted by specified OUTRIDER groups
    """
    def createGroupIds(self, group_key="DROP_GROUP", assay_key="RNA_ASSAY", sep=','):
        
        # Get unique groups
        ids = self.getSampleIDs(assay_key)
        df = self.sample_annotation[self.sample_annotation[assay_key].isin(ids)]
        df = df[[assay_key, group_key]].drop_duplicates().copy()
        
        # get unique group names
        groups = []
        for s in set(df[group_key]):
            groups.extend(s.split(sep))
        groups = set(groups)
        
        # collect IDs per group
        return {gr : df[df[group_key].str.contains(f'(^|{sep}){gr}({sep}|$)')][assay_key].tolist() 
                        for gr in groups}
    
    def subsetGroups(self, ids_by_group, subset_groups):
        if subset_groups is None:
            return ids_by_group
        else:
            return {gr:ids for gr, ids in ids_by_group.items() if gr in subset_groups}
    
    def getGeneAnnotationFile(self, annotation):
        return self.config["GENE_ANNOTATION"][annotation]

 
        

