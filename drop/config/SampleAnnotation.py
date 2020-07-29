import pandas as pd
from pathlib import Path
from snakemake.logging import logger
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')

class SampleAnnotation:
    
    SAMPLE_ANNOTATION_COLUMNS = ["RNA_ID", "RNA_BAM_FILE", "DNA_ID", "DNA_VCF_FILE",
                                 "DROP_GROUP","PAIRED_END", "COUNT_MODE",
                                 "COUNT_OVERLAPS", "STRAND"]
    
    def __init__(self, file, root):
        """
        sa_file: sample annotation file location from config
        root: output location for file mapping
        """
        self.root = Path(root)
        self.file = file
        self.sa = self.parse()
        self.idMapping = self.createIdMapping()
        self.sampleFileMapping = self.createSampleFileMapping()
        
        self.rnaIDs = self.createGroupIds(group_key="DROP_GROUP",
                                          file_type="RNA_BAM_FILE", sep=',')
        self.dnaIDs = self.createGroupIds(group_key="DROP_GROUP",
                                          file_type="DNA_VCF_FILE", sep=',')
        
    def parse(self, sep='\t'):
        """
        read and check sample annotation for missing columns
        clean columns and set types
        """
        sa = pd.read_csv(self.file, sep=sep)
        missing_cols = [x for x in self.SAMPLE_ANNOTATION_COLUMNS if x not in sa.columns.values]
        if len(missing_cols) > 0:
            raise ValueError("Incorrect columns in sample annotation. Missing:"+
                             f"\n {missing_cols}\n")
        
        # remove unwanted characters
        col = sa["DROP_GROUP"]
        col = col.str.replace(" ", "").str.replace("(|)", "", regex=True)
        sa["DROP_GROUP"] = col
        
        # set ID type as string
        for ID_key in ["RNA_ID", "DNA_ID"]:
            sa[ID_key] = sa[ID_key].apply(
                    lambda x: str(x) if not pd.isnull(x) else x
            )
        return sa
    
    def createIdMapping(self):
        """
        Get mapping of RNA and DNA IDs
        """
        return self.sa[["RNA_ID", "DNA_ID"]].drop_duplicates().dropna()
    
    def createSampleFileMapping(self):
        """
        create a sample file mapping with unique entries of existing files
            columns: [ID | ASSAY | FILE_TYPE | FILE_PATH ]
        """

        assay_mapping = {'RNA_ID': ['RNA_BAM_FILE', 'GENE_COUNTS_FILE'], 'DNA_ID': ['DNA_VCF_FILE']}
        assay_subsets = []
        for id_, file_types in assay_mapping.items():
            for file_type in file_types:
                df = self.sa[[id_, file_type]].dropna().drop_duplicates().copy()
                df.rename(columns={id_: 'ID', file_type: 'FILE_PATH'}, inplace=True)
                df['ASSAY'] = id_
                df['FILE_TYPE'] = file_type
                assay_subsets.append(df)
        file_mapping = pd.concat(assay_subsets)
        
        # cleaning SAMPLE_FILE_MAPPING
        file_mapping.dropna(inplace=True)
        file_mapping.drop_duplicates(inplace = True)
        
        # check for missing files
        exist = [Path(x).exists() for x in file_mapping["FILE_PATH"]]
        if sum(exist) == 0:
            message = "File mapping is empty. "
            message += "Please check that all files in your sample annotation exist."
            raise ValueError(message)
        elif sum(exist) < file_mapping.shape[0]:
            file_mapping = file_mapping[exist]
            missing = [x for x in file_mapping["FILE_PATH"] if not Path(x).exists()]
            info = f"WARNING: {len(missing)} "
            info += "of the files in the samples annotation do not exist"
            logger.info(info)
            logger.debug(f"missing files: {missing}")
        
        file_mapping.to_csv(self.root / "file_mapping.csv", index=False)

        return file_mapping
    
    def subsetFileMapping(self, file_type = None, sample_id = None):
        """
        subset by one or more values of different columns from sample file mapping
            file_type: file type/types, corresponding to 'FILE_TYPE' column
            sample_id: sample ID/IDs
        """
        
        def subsetBy(df, values, column):
            if values is None:
                return df
            elif isinstance(values, str):
                return df[df[column] == values]
            else:
                return df[df[column].isin(values)]
        
        subset = self.sampleFileMapping
        subset = subsetBy(subset, file_type, "FILE_TYPE")
        subset = subsetBy(subset, sample_id, "ID")
        
        return subset
    
    def getFilePath(self, sample_id, file_type, single_file=True):
        """
        Get path to input data file by sample ID
        """
        path = self.subsetFileMapping(file_type, sample_id)["FILE_PATH"]
        path = path.tolist()
        if single_file:
            if len(path) > 1:
                message = "Trying to return more than 1 path for 1 sample ID and file type"
                raise ValueError(message)
            path = path[0]
        return path
    
    def getFilePaths(self, file_type, group=None):
        """
        Get all file paths of a file type
            file_type: 'RNA_BAM_FILE' or 'DNA_VCF_FILE'
            group: name of DROP_GROUP
        """
        if group is None:
            sampleIDs = self.getSampleIDs(file_type)
        else:
            sampleIDs = self.getGroupedIDs(file_type)[group]
        return self.getFilePath(sampleIDs, file_type, single_file=False)
    
    ### DROP Groups ###
    
    def getGroupedIDs(self, assay):
        """
        assay: what assay the IDs should be from. Can be file_type or 'RNA'/'DNA'
        """
        if "RNA" in assay:
            return self.rnaIDs
        elif "DNA" in assay:
            return self.dnaIDs
        
    def getGroups(self, assay="RNA"):
        return self.getGroupedIDs(assay).keys()
    
    def getIDsByGroup(self, group, assay="RNA"):
        return self.getGroupedIDs(assay)[group]
    
    def createGroupIds(self, group_key="DROP_GROUP", file_type="RNA_BAM_FILE", sep=','):
        """
        Create a mapping of DROP groups to sample IDs
        """
        sa = self.sa
        sf = self.sampleFileMapping

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
        grouped = {gr : df[df[group_key].str.contains(f'(^|{sep}){gr}({sep}|$)')][assay_id].tolist() 
                        for gr in groups}
        # remove groups labeled as None
        grouped = {gr : list(set(ids)) for gr, ids in grouped.items() if gr is not None}
        return grouped
    
    def subsetGroups(self, subset_groups, assay="RNA", warn=30, error=10):
        """
        warn : number of samples threshold at which to warn about too few samples
        error: number of samples threshold at which to give error
        """
        ids_by_group = self.getGroupedIDs(assay)
        
        if subset_groups is None:
            subset = ids_by_group
        else:
            subset_groups = [subset_groups] if subset_groups.__class__ == str else subset_groups
            subset = {gr:ids for gr, ids in 
                      ids_by_group.items() if gr in subset_groups}
        
        for group in subset_groups:
            if len(subset[group]) < error:
                message = f'Too few IDs in DROP_GROUP {group}'
                message += f', please ensure that it has at least {error} IDs'
                message += f', groups: {subset[group]}'
                raise ValueError(message)
            elif len(subset[group]) < warn:
                logger.info(f'WARNING: Less than {warn} IDs in DROP_GROUP {group}')
        
        return subset
    
    def getSampleIDs(self, file_type):
        ids = self.subsetFileMapping(file_type)["ID"]
        return list(ids)
    
    def checkFileExists(self, sampleID, file_type, verbose=True):
        """
        note: we already checked for non-existing files in the init, so we
        only need to check whether the ID is in the sample_file_mapping here
        """ 
        x = self.sampleFileMapping.query("(FILE_TYPE == @file_type) & (ID == @sampleID)")["FILE"]
        exists = (len(x) != 0)
        if (not exists) and verbose:
            logger.debug(f"FILE NOT FOUND FOR sampleID: {sampleID} and file type {file_type}")
        return exists
    