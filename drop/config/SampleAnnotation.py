from drop import utils
import pandas as pd
from pathlib import Path
from collections import defaultdict
from snakemake.logging import logger
import warnings

warnings.filterwarnings("ignore", 'This pattern has match groups')


class SampleAnnotation:
    FILE_TYPES = ["RNA_BAM_FILE", "DNA_VCF_FILE", "GENE_COUNTS_FILE"]
    SAMPLE_ANNOTATION_COLUMNS = FILE_TYPES + [
        "RNA_ID", "DNA_ID", "DROP_GROUP", "ANNOTATION",
        "PAIRED_END", "COUNT_MODE", "COUNT_OVERLAPS", "STRAND"
    ]

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

        self.rnaIDs = self.createGroupIds(file_type="RNA_BAM_FILE", sep=',')
        self.dnaIDs = self.createGroupIds(file_type="DNA_VCF_FILE", sep=',')
        # external counts
        self.extGeneCountIDs = self.createGroupIds(file_type="GENE_COUNTS_FILE", sep=',')

    def parse(self, sep='\t'):
        """
        read and check sample annotation for missing columns
        clean columns and set types
        """
        data_types = {
            "RNA_ID": str, "DNA_ID": str, "DROP_GROUP": str, "ANNOTATION": str,
            "PAIRED_END": bool, "COUNT_MODE": str, "COUNT_OVERLAPS": bool, "STRAND": str
        }
        sa = pd.read_csv(self.file, sep=sep, index_col=False, converters=data_types)
        missing_cols = [x for x in self.SAMPLE_ANNOTATION_COLUMNS if x not in sa.columns.values]
        if len(missing_cols) > 0:
            raise ValueError(f"Incorrect columns in sample annotation file. Missing:\n{missing_cols}")

        # remove unwanted characters
        sa["DROP_GROUP"] = sa["DROP_GROUP"].str.replace(" ", "").str.replace("(|)", "", regex=True)

        return sa

    ##### Construction

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
        file_mapping.drop_duplicates(inplace=True)

        # check for missing files
        existing = utils.checkFileExists(file_mapping["FILE_PATH"])
        if len(existing) == 0:
            message = "File mapping is empty. "
            message += "Please check that all files in your sample annotation exist."
            raise FileNotFoundError(message)
        elif len(existing) < file_mapping.shape[0]:
            missing = set(file_mapping["FILE_PATH"]) - set(existing)
            logger.info(f"WARNING: {len(missing)} files missing in samples annotation. Ignoring...")
            logger.debug(f"Missing files: {missing}")
            file_mapping = file_mapping[file_mapping["FILE_PATH"].isin(existing)]

        # write file mapping
        file_mapping.to_csv(self.root / "file_mapping.csv", index=False)
        return file_mapping

    def createGroupIds(self, group_key="DROP_GROUP", file_type=None, sep=','):
        """
        :param group_key: name of group column in sample annotation
        :param file_type: name of file column e.g. "RNA_BAM_FILE", "DNA_VCF_FILE"
        :param sep: separator of multiple groups in group column
        :return: mapping of drop group and ID
        """
        if not file_type:
            file_type = "RNA_BAM_FILE"
        # infer ID type from file type
        assay_id = self.subsetFileMapping(file_type)["ASSAY"].unique()
        if len(assay_id) == 0:
            return defaultdict(list)  # no files, return empty mapping
        elif len(assay_id) > 1:
            raise ValueError(f"More than 1 assay entry for file type {file_type}:\n{assay_id}")
        else:
            assay_id = assay_id[0]

        # Subset sample annotation to only IDs of specified file_type
        ids = self.getSampleIDs(file_type)
        df = self.sa[self.sa[assay_id].isin(ids)]
        # mapping of ID to group names
        df = df[[assay_id, group_key]].drop_duplicates().copy()

        # get unique group names
        groups = []
        for s in set(self.sa[group_key]):
            groups.extend(s.split(sep))
        groups = set(groups)

        # collect IDs per group
        grouped = {gr: df[df[group_key].str.contains(f'(^|{sep}){gr}({sep}|$)')][assay_id].tolist()
                   for gr in groups}
        # remove groups labeled as None
        grouped = {gr: list(set(ids)) for gr, ids in grouped.items() if gr is not None}
        return grouped

    ### Subsetting

    def subsetSampleAnnotation(self, column, values, subset=None):
        """
        subset by one or more values of different columns from sample file mapping
            :param column: valid column in sample annotation
            :param values: values of column to subset
            :param subset: subset sample annotation
        """
        sa_cols = set(self.SAMPLE_ANNOTATION_COLUMNS)
        if subset is None:
            subset = self.sa
        else:
            # check type for subset
            if not isinstance(subset, pd.DataFrame):
                raise TypeError(f"Is not pandas DataFrame\n {subset}")
            if not sa_cols <= set(subset.columns):  # check if mandatory cols not contained
                raise ValueError(f"Subset columns not the same as {sa_cols}\ngot: {subset.columns}")

        # check if column is valid
        if not column in sa_cols:
            raise KeyError(f"Column '{column}' invalid for sample annotation")
        return utils.subsetBy(subset, column, values)

    def subsetFileMapping(self, file_type=None, sample_id=None):
        """
        subset by one or more values of different columns from sample file mapping
            file_type: file type/types, corresponding to 'FILE_TYPE' column
            sample_id: sample ID/IDs
        """
        subset = self.sampleFileMapping
        subset = utils.subsetBy(subset, "FILE_TYPE", file_type)
        subset = utils.subsetBy(subset, "ID", sample_id)
        return subset

    def subsetGroups(self, subset_groups, assay="RNA"):
        """
        Subset DROP group to sample IDs mapping by list of groups (`subset_groups`).
        :param subset_groups: list of groups to include
        :param assay: name/prefix of assay type
        :return: dictionary with group names as keys and ID lists as entries
        """
        ids_by_group = self.getGroupedIDs(assay)

        if subset_groups is None:
            subset = ids_by_group
        else:
            subset_groups = [subset_groups] if subset_groups.__class__ == str else subset_groups
            subset = {gr: ids for gr, ids in
                      ids_by_group.items() if gr in subset_groups}
            subset = defaultdict(list, subset)
        return subset

    ### Getters

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

    def getImportCountFiles(self, annotation, group, file_type="GENE_COUNTS_FILE",
                            annotation_key="ANNOTATION", group_key="DROP_GROUP"):
        """
        :param annotation: annotation name as specified in config and ANNOTATION column
        :param group: a group of the DROP_GROUP column
        :return: set of unique external count file names
        """
        subset = self.subsetSampleAnnotation(annotation_key, annotation)
        subset = self.subsetSampleAnnotation(group_key, group, subset)
        return set(subset[file_type].tolist())

    def getRow(self, column, value):
        sa = self.sa
        if column not in sa.columns:
            raise KeyError(f"column {column} not in sample annotation")
        row = sa[sa[column] == value]
        if row.shape[0] != 1:
            raise ValueError(f"sa[sa[{column}] == {value}] should have 1 row")
        return row

    ### DROP Groups ###

    def getGroupedIDs(self, assays):
        """
        Get group to IDs mapping
        :param assays: list of or single assay the IDs should be from. Can be file_type or 'RNA'/'DNA'
        """
        assays = [assays] if isinstance(assays, str) else assays
        groupedIDs = defaultdict(list)
        for assay in assays:
            if "RNA" in assay:
                groupedIDs.update(self.rnaIDs)
            elif "DNA" in assay:
                groupedIDs.update(self.dnaIDs)
            elif "GENE_COUNT" in assay:
                groupedIDs.update(self.extGeneCountIDs)
            else:
                raise ValueError(f"'{assay}' is not a valid assay name")
        return groupedIDs

    def getGroups(self, assay="RNA"):
        return self.getGroupedIDs(assay).keys()

    def getIDsByGroup(self, group, assay="RNA"):
        try:
            ids = self.getGroupedIDs(assay)[group]
        except:
            ids = []
        return ids

    def getSampleIDs(self, file_type):
        ids = self.subsetFileMapping(file_type)["ID"]
        return list(ids)
