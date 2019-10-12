import os
import pandas as pd
import wbuild
import pathlib
from snakemake.logging import logger
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')

class ConfigHelper:
    
    def __init__(self, config):
        
        self.config = self.setDefaults(config)
        
        # sample annotation
        self.sample_annotation = self.getSampleAnnotation(config["sampleAnnotation"])
        self.sample_file_mapping = self.createSampleFileMapping(self.sample_annotation)
        
        # Group IDs
        # remove unwanted characters
        self.sample_annotation["DROP_GROUP"] = self.sample_annotation["DROP_GROUP"].str.replace("(", "").str.replace(")", "")
        self.all_rna_ids = self.createGroupIds(group_key="DROP_GROUP", file_type="RNA_BAM_FILE", sep=',')
        
        ## outrider
        self.config["useGeneNames"] = self.config["aberrantExpression"]["useGeneNames"]
        
        if self.config["aberrantExpression"]["groups"] is None:
            self.config["aberrantExpression"]["groups"] = list(self.all_rna_ids.keys())

        if self.config["aberrantExpression"]["minIds"] is None:
            self.config["aberrantExpression"]["minIds"] = 40
        
        self.outrider_all = self.subsetGroups(self.all_rna_ids, self.config["aberrantExpression"]["groups"])
        self.outrider_filtered = {name:ids 
            for name, ids in self.outrider_all.items()
            if len(ids) > self.config["aberrantExpression"]["minIds"]
        }
        self.config["outrider_all"], self.config["outrider_filtered"] = self.outrider_all, self.outrider_filtered
        
        ## fraser
        if self.config["aberrantSplicing"]["groups"] is None:
            self.config["aberrantSplicing"]["groups"] = list(self.all_rna_ids.keys())
        
        if self.config["aberrantSplicing"]["minIds"] is None:
            self.config["aberrantSplicing"]["minIds"] = 40
        self.fraser_all = self.subsetGroups(self.all_rna_ids, self.config["aberrantSplicing"]["groups"])
        self.fraser_filtered = {name:ids 
            for name, ids in self.fraser_all.items()
            if len(ids) > self.config["aberrantSplicing"]["minIds"]
        }
        self.config["fraser_all"], self.config["fraser_filtered"] = self.fraser_all, self.fraser_filtered
        
        ## mae
        if self.config["mae"]["groups"] is None:
            self.config["mae"]["groups"] = list(self.all_rna_ids.keys())
        mae_rna_by_group = self.subsetGroups(self.all_rna_ids, self.config["mae"]["groups"])
        self.mae_ids = self.createMaeIDS(mae_rna_by_group, id_sep='--')
        self.config["mae_ids"] = self.mae_ids
        
    def setDefaults(self, config):
        """
        set defaults for config keys
        """
        
        if config is None:
            wconf = wbuild.utils.Config()
            config = wconf.conf_dict
        
        # check for missing keys
        def check_keys(dict_, keys):
            for key in keys:
                try:
                    dict_[key]
                except:
                    raise KeyError(f"{key} is mandatory but missing")
                    
        mandatory_keys = ["htmlOutputPath", "root", "geneAnnotation", "sampleAnnotation",
                           "aberrantExpression", "aberrantSplicing", "mae"]
        check_keys(config, keys=mandatory_keys)
        check_keys(config["mae"], keys=["genome", "qcVcf"])
        check_keys(config["mae"]["qcVcf"], keys=["UCSC", "NCBI"])
        
        # set default
        def setKey(dict_, sub, key, default):
            if sub is not None:
                for x in sub:
                    dict_ = dict_[x]
            if key not in dict_ or dict_[key] is None:
                print(f'{key} not in config{sub}, using default')
                dict_[key] = default
        
        config["indexWithFolderName"] = True
        config["fileRegex"] = ".*\.R"
        setKey(config, None, "projectTitle", "DROP: Detection of RNA Outlier Pipeline")
        setKey(config, None, "tmpdir", os.path.join(config["root"], "tmp"))
        setKey(config, ["aberrantExpression"], "fpkmCutoff", 1)
        setKey(config, ["aberrantExpression"], "groups", None)
        setKey(config, ["aberrantExpression"], "padjCutoff", .05)
        setKey(config, ["aberrantExpression"], "zscoreCutoff", 0)
        setKey(config, ["aberrantExpression"], "useGeneNames", True)
        setKey(config, ["aberrantSplicing"], "groups", None)
        setKey(config, ["aberrantSplicing"], "deltaPsiCutoff", 0.05)
        setKey(config, ["mae"], "geneAssembly", "hg19")
        setKey(config, ["mae"], "gatkIgnoreHeaderCheck", True)
        setKey(config, ["mae"], "padjCutoff", .05)
        setKey(config, ["mae"], "allelicRatioCutoff", 0.8)
        setKey(config, ["mae"], "maxAF", .001)
        setKey(config, ["mae"], "groups", None)
        setKey(config, ["mae"], "qcGroup", None)
        setKey(config, ["tools"], "samtoolsCmd", "samtools")
        setKey(config, ["tools"], "bcftoolsCmd", "bcftools")
        setKey(config, ["tools"], "gatkCmd", "gatk")
        
        return config
        
    def getSampleAnnotation(self, filename, sep='\t'):
        """
        read and check sample annotation for missing columns
        """
        sample_annotation = pd.read_csv(filename, sep=sep)
        #TODO add more column names
        mandatory_columns = {"RNA_ID", "RNA_BAM_FILE", "DNA_ID", "DNA_VCF_FILE"}
        missing_cols = [x for x in mandatory_columns if x not in sample_annotation.columns.values]
        if len(missing_cols) > 0:
            raise ValueError(f"Incorrect columns in sample annotation, missing:\n {missing_cols}\n")
        return sample_annotation
    
    def createSampleFileMapping(self, sample_annotation):
        """
        create a sample file mapping with unique entries of existing files
            columns: [ID | ASSAY | FILE_TYPE | FILE_PATH ]
        """

        assay_mapping = {'RNA_ID':'RNA_BAM_FILE', 'DNA_ID':'DNA_VCF_FILE'}                        
        
        assay_subsets = []
        for id_, file_type in assay_mapping.items():
            df = sample_annotation[[id_, file_type]].drop_duplicates().copy()
            df.rename(columns={id_:'ID', file_type:'FILE_PATH'}, inplace=True)
            df['ASSAY'] = id_
            df['FILE_TYPE'] = file_type
            assay_subsets.append(df)
        
        file_mapping = pd.concat(assay_subsets)
        
        # cleaning SAMPLE_FILE_MAPPING
        file_mapping.dropna(inplace=True)
        existent = [os.path.exists(x) for x in file_mapping["FILE_PATH"]]
        file_mapping = file_mapping[existent].drop_duplicates()
        
        file_mapping.to_csv(self.getProcDataDir() + "/file_mapping.csv", index=False)
        
        return file_mapping
        
    """ 
    Get directory path for processed data
    """
    def getProcDataDir(self):
        return self.config["root"] + "/processed_data"
    
    """ 
    Get directory path for processed results
    """
    def getProcResultsDir(self):
        return self.config["root"] + "/processed_results"
    
    """
    Get sample ID by file type
    """
    def getSampleIDs(self, file_type):
        fm = self.sample_file_mapping
        return list(fm[fm["FILE_TYPE"] == file_type]["ID"]) 
    
    
    def checkFileExists(self, sampleID, file_type, verbose=True):
        # note: we already checked for non-existing files in the init, so we
        # only need to check whether the ID is in the sample_file_mapping here
        x = self.sample_file_mapping.query("(FILE_TYPE == @file_type) & (ID == @sampleID)")["FILE"]
        exists = (len(x) != 0)
        if (not exists) and verbose:
            print(f"FILE NOT FOUND FOR sampleID: {sampleID} and file type {file_type}")
        return exists
        
    """
    Returns vcf and rna files for MAE pipeline
    """
    def getMaeByGroup(self, group):
        if not isinstance(group, str):
            group = list(group)[0]
        return self.mae_ids[group]

    def getMaeAll(self):
        all_ids = []
        for group in self.config["mae"]["groups"]:
            all_ids.extend(self.mae_ids[group])
        return all_ids

    def createMaeIDS(self, rna_id_by_group, id_sep='--'):
        
        all_mae_files = self.getAllMaeFiles()
        
        # subset by group
        mae_ids = {}
        for gr, rna_ids in rna_id_by_group.items():
            mae_subset = all_mae_files [all_mae_files ["RNA_ID"].isin(rna_ids)]
            vcf_rna_pairs = zip(mae_subset["DNA_ID"], mae_subset["RNA_ID"])
            mae_ids[gr] = list(map(id_sep.join, vcf_rna_pairs))
        
        return mae_ids
        
    def getAllMaeFiles(self):
        
        mae_files = self.sample_annotation[["RNA_ID", "DNA_ID"]]
        mae_files = mae_files.dropna()
        
        return mae_files
    
    def getRNAByGroup(self, group):
        if not isinstance(group, str):
            group = list(group)[0]
        return self.all_rna_ids[group]
    
    def getFilePath(self, sampleId, file_type):
        """
        Function for getting the file path given the sampleId and file t
        sampleId: ID of sample
        file_type: e.g. "RNA_BAM_FILE", "DNA_VCF_FILE"
        """
        fm = self.sample_file_mapping
        if isinstance(file_type, str):
            path = fm.query("FILE_TYPE == @file_type")
        else:
            path = fm[fm["FILE_TYPE"].isin(file_type)]
        path = path.query("ID == @sampleId")["FILE_PATH"]
        return path.iloc[0]
    
    def getFilePaths(self, file_type, group=None, ids_by_group=None):
        """
        file_type: e.g. "RNA_BAM_FILE", "DNA_VCF_FILE"
        group: name of dataset/ group
        ids_by_group: dictionary of IDs by group names
        """
        
        # subset by group if group is specified
        if group is None or ids_by_group is None:
            sampleIDs = self.sample_file_mapping.query("FILE_TYPE == @file_type")["ID"]
        else:
            sampleIDs = ids_by_group[group]
        
        files = [] # collect file names
        for sampleID in sampleIDs:
            files.append(self.getFilePath(sampleID, file_type))
        return files
    
    def createGroupIds(self, group_key="DROP_GROUP", file_type="RNA_BAM_FILE", sep=','):
        """
        Create a full and filtered list of RNA assay IDs subsetted by specified OUTRIDER groups
        """
        sa = self.sample_annotation
        sf = self.sample_file_mapping
        
        assay_id = sf[sf["FILE_TYPE"] == file_type]["ASSAY"].iloc[0]
        
        # Get unique groups
        ids = self.getSampleIDs(file_type)
        df = sa[sa[assay_id].isin(ids)]
        df = df[[assay_id, group_key]].drop_duplicates().copy()
        
        # get unique group names
        groups = []
        for s in set(sa[group_key]):
            groups.extend(s.split(sep))
        groups = set(groups)
        
        # collect IDs per group
        return {gr : df[df[group_key].str.contains(f'(^|{sep}){gr}({sep}|$)')][assay_id].tolist() 
                        for gr in groups}
    
    def subsetGroups(self, ids_by_group, subset_groups):
        if subset_groups is None:
            return ids_by_group
        else:
            return {gr:ids for gr, ids in ids_by_group.items() if gr in subset_groups}
    
    def getGeneAnnotationFile(self, annotation):
        return self.config["geneAnnotation"][annotation]

 
        

