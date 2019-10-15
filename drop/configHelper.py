import os
import pandas as pd
import wbuild
import pathlib
from snakemake.logging import logger
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')

#TODO add more column names
SAMPLE_ANNOTATION_COLUMNS = ["RNA_ID", "RNA_BAM_FILE", "DNA_ID", "DNA_VCF_FILE", "DROP_GROUP",
"PAIRED_END", "COUNT_MODE", "COUNT_OVERLAPS", "STRAND"]

class ConfigHelper:
    
    def __init__(self, config):
        
        if config is None:
            wconf = wbuild.utils.Config()
            config = wconf.conf_dict
        
        self.config = self.setDefaults(config)
        self.createDirs()
        
        # sample annotation
        self.sample_annotation = self.getSampleAnnotation(config["sampleAnnotation"])
        self.sample_file_mapping = self.createSampleFileMapping(self.sample_annotation)

        self.all_rna_ids = self.createGroupIds(group_key="DROP_GROUP", file_type="RNA_BAM_FILE", sep=',')
        
        # aberrantExpression
        groups = self.setKey(self.config, ["aberrantExpression"], "groups", self.all_rna_ids.keys())
        self.outrider_ids = self.subsetGroups(self.all_rna_ids, groups)
        self.config["outrider_ids"] = self.outrider_ids
        
        # aberrantSplicing
        groups = self.setKey(self.config, ["aberrantSplicing"], "groups", self.all_rna_ids.keys())
        self.fraser_ids = self.subsetGroups(self.all_rna_ids, groups)
        self.config["fraser_ids"] = self.fraser_ids
        
        # mae
        groups = self.setKey(self.config, ["mae"], "groups", self.all_rna_ids.keys())
        self.mae_ids = self.createMaeIDS(self.all_rna_ids, groups, id_sep='--')
        self.config["mae_ids"] = self.mae_ids
        
        
    def createDirs(self):
        
        def createIfMissing(directory):
            if not os.path.exists(directory):
                print(f"creating {directory}")
                os.makedirs(directory)
        
        createIfMissing(self.getProcDataDir())
        createIfMissing(self.getProcResultsDir())
        
    def checkConfig(self, config):
        # check for missing keys
        def check_keys(dict_, keys):
            for key in keys:
                try:
                    dict_[key]
                except:
                    raise KeyError(f"{key} is mandatory but missing")
                    
        FILE_KEYS = ["htmlOutputPath", "root", "geneAnnotation", "sampleAnnotation", "mae"]
        check_keys(config, keys=FILE_KEYS)
        check_keys(config["mae"], keys=["genome", "qcVcf"])
        check_keys(config["mae"]["qcVcf"], keys=["UCSC", "NCBI"])
    
    def setDefaults(self, config):
        """
        set defaults for config keys
        """
        
        config["indexWithFolderName"] = True
        config["fileRegex"] = ".*\.R"
        
        setKey = self.setKey
        setKey(config, None, "projectTitle", "DROP: Detection of RNA Outlier Pipeline")
        setKey(config, None, "tmpdir", os.path.join(config["root"], "tmp"))
        
        # aberrant expression
        setKey(config, None, "aberrantExpression", dict())
        setKey(config, ["aberrantExpression"], "fpkmCutoff", 1)
        setKey(config, ["aberrantExpression"], "groups", None)
        setKey(config, ["aberrantExpression"], "padjCutoff", .05)
        setKey(config, ["aberrantExpression"], "zscoreCutoff", 0)
        setKey(config, ["aberrantExpression"], "useGeneNames", True)
        
        # aberrant splicing
        setKey(config, None, "aberrantSplicing", dict())
        setKey(config, ["aberrantSplicing"], "groups", None)
        setKey(config, ["aberrantSplicing"], "deltaPsiCutoff", 0.05)
        
        # monoallelic expression
        setKey(config, None, "mae", dict())
        setKey(config, ["mae"], "geneAssembly", "hg19")
        setKey(config, ["mae"], "gatkIgnoreHeaderCheck", True)
        setKey(config, ["mae"], "padjCutoff", .05)
        setKey(config, ["mae"], "allelicRatioCutoff", 0.8)
        setKey(config, ["mae"], "maxAF", .001)
        setKey(config, ["mae"], "groups", None)
        setKey(config, ["mae"], "qcGroup", None)
        
        # commandline tools
        setKey(config, None, "tools", dict())
        setKey(config, ["tools"], "samtoolsCmd", "samtools")
        setKey(config, ["tools"], "bcftoolsCmd", "bcftools")
        setKey(config, ["tools"], "gatkCmd", "gatk")
        
        return config
    
    def setKey(self, dict_, sub, key, default):
        if sub is not None:
            for x in sub:
                dict_ = dict_[x]
        if key not in dict_ or dict_[key] is None:
            print(f'{key} not in config{sub}, using default')
            dict_[key] = default
        return dict_[key]
    
    def getSampleAnnotation(self, filename, sep='\t'):
        """
        read and check sample annotation for missing columns
        """
        sample_annotation = pd.read_csv(filename, sep=sep)
        missing_cols = [x for x in SAMPLE_ANNOTATION_COLUMNS if x not in sample_annotation.columns.values]
        if len(missing_cols) > 0:
            raise ValueError(f"Incorrect columns in sample annotation, missing:\n {missing_cols}\n")
            
        # remove unwanted characters
        col = sample_annotation["DROP_GROUP"]
        col = col.str.replace(" ", "").str.replace("(|)", "", regex=True)
        sample_annotation["DROP_GROUP"] = col
        
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
        if sum(existent) < file_mapping.shape[0]:
            print("WARNING: there are files in the sample annotation that do not exist")
        file_mapping = file_mapping[existent].drop_duplicates()
        if file_mapping.shape[0] == 0:
            raise ValueError("No files exist in sample annotation. Please check your sample annotation.")
        
        file_mapping.to_csv(self.getProcDataDir() + "/file_mapping.csv", index=False)
        
        return file_mapping
    
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
    
    def subsetGroups(self, ids_by_group, subset_groups, warn=30, error=10):
        if subset_groups is None:
            subset = ids_by_group
        else:
            subset = {gr:ids for gr, ids in ids_by_group.items() if gr in subset_groups}
        
        for group in subset_groups:
            if len(subset[group]) < error:
                raise ValueError(f'Too few IDs in DROP_GROUP {group}, please ensure that it has at least {error} IDs')
            elif len(subset[group]) < warn:
                print(f'WARNING: Less than {warn} IDs in DROP_GROUP {group}')
        
        return subset
    
    def createMaeIDS(self, ids_by_group, subset_groups, id_sep='--'):
        """
        create MAE IDs from smaple annotation
        """
        rna_id_by_group = self.subsetGroups(ids_by_group, subset_groups)
        all_mae_files = self.sample_annotation[["RNA_ID", "DNA_ID"]].drop_duplicates().dropna()
        
        # subset by group
        mae_ids = {}
        for gr, rna_ids in rna_id_by_group.items():
            mae_subset = all_mae_files [all_mae_files ["RNA_ID"].isin(rna_ids)]
            vcf_rna_pairs = zip(mae_subset["DNA_ID"], mae_subset["RNA_ID"])
            mae_ids[gr] = list(map(id_sep.join, vcf_rna_pairs))
        
        return mae_ids
    
    
    ## GETTERS
    def getProcDataDir(self):
        return self.config["root"] + "/processed_data"
    
    def getProcResultsDir(self):
        return self.config["root"] + "/processed_results"
    
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
        group: name of DROP_GROUP
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

    def getRNAByGroup(self, group):
        if not isinstance(group, str):
            group = list(group)[0]
        return self.all_rna_ids[group]
    
    def getMaeByGroup(self, group):
        if not isinstance(group, str):
            group = list(group)[0]
        return self.mae_ids[group]

    def getMaeAll(self):
        """
        Get a list of all MAE IDs from the groups specified in the config.
        Useful for collecting all MAE IDs ungrouped.
        """
        all_ids = []
        for group in self.config["mae"]["groups"]:
            all_ids.extend(self.mae_ids[group])
        return all_ids
    
    def getGeneAnnotationFile(self, annotation):
        return self.config["geneAnnotation"][annotation]

 
        

