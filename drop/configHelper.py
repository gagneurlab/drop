from drop import submodules
import pandas as pd
import wbuild
import pathlib
from snakemake.logging import logger
from snakemake.io import expand
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')

SAMPLE_ANNOTATION_COLUMNS = ["RNA_ID", "RNA_BAM_FILE", "DNA_ID", "DNA_VCF_FILE", "DROP_GROUP",
"PAIRED_END", "COUNT_MODE", "COUNT_OVERLAPS", "STRAND"]
VERBOSE = True

class ConfigHelper:
    
    def __init__(self, config, method=None):
        
        self.method = method
        if config is None:
            wconf = wbuild.utils.Config()
            config = wconf.conf_dict
        
        config = self.checkConfig(config)
        self.config = self.setDefaults(config)
        self.createDirs()
    
    # Function to remove duplicate elements inside a DROP group    
    def drop_dups(self, l):
      return(list(set(l)))
      
    def parse(self):
        """
        parse sample annotation, create sample-file mapping and extract IDs for each submodule
        """
        # sample annotation
        self.sample_annotation = self.getSampleAnnotation(self.config["sampleAnnotation"])
        self.sample_file_mapping = self.createSampleFileMapping(self.sample_annotation)

        self.all_rna_ids = self.createGroupIds(group_key="DROP_GROUP", file_type="RNA_BAM_FILE", sep=',')
        
        if self.method == "AE" or self.method is None:
            groups = self.setKey(self.config, ["aberrantExpression"], "groups", self.all_rna_ids.keys())
            self.outrider_ids = self.subsetGroups(self.all_rna_ids, groups)
            for k,v in self.outrider_ids.items():
              self.outrider_ids[k] = self.drop_dups(v)
            self.config["outrider_ids"] = self.outrider_ids
        
        if self.method == "AS" or self.method is None:
            groups = self.setKey(self.config, ["aberrantSplicing"], "groups", self.all_rna_ids.keys())
            self.fraser_ids = self.subsetGroups(self.all_rna_ids, groups)
            for k,v in self.fraser_ids.items():
              self.fraser_ids[k] = self.drop_dups(v)
            self.config["fraser_ids"] = self.fraser_ids
        
        if self.method == "MAE" or self.method is None:
            groups = self.setKey(self.config, ["mae"], "groups", self.all_rna_ids.keys())
            self.config["mae"]["qcGroups"] = self.setKey(self.config, ["mae"], "qcGroups", groups)
            self.mae_ids = self.createMaeIDS(self.all_rna_ids, groups, id_sep='--')
            for k,v in self.mae_ids.items():
              self.mae_ids[k] = self.drop_dups(v)
            self.config["mae_ids"] = self.mae_ids
        
        return self.config

    
    #### FILE CHECKING AND DEFAULT SETTINGS
    
    def createDirs(self):
        
        def createIfMissing(directory):
            directory = pathlib.Path(directory)
            if not directory.exists():
                logger.debug(f"creating {directory}")
                directory.mkdir(parents=True)
        
        createIfMissing(self.getProcDataDir())
        createIfMissing(self.getProcResultsDir())
        
    def checkConfig(self, config):
        
        def check_keys(dict_, keys):
            keys = dict_.keys() if keys is None else keys
            for key in keys:
                try:
                    dict_[key]
                except:
                    raise KeyError(f"{key} is mandatory but missing")
                # get real path
                if isinstance(dict_[key], str):
                    filename = dict_[key]
                    dict_[key] = str(pathlib.Path(filename).expanduser())
        
        FILE_KEYS = ["htmlOutputPath", "root", "geneAnnotation", "sampleAnnotation", "mae"]
        check_keys(config, keys=FILE_KEYS)
        check_keys(config["geneAnnotation"], keys=None)
        check_keys(config["mae"], keys=["genome", "qcVcf"])
        
        return config
    
    def setDefaults(self, config, method=None):
        """
        set defaults for config keys
        """
        
        config["indexWithFolderName"] = True
        config["fileRegex"] = ".*\.R"
        
        setKey = self.setKey
        setKey(config, None, "projectTitle", "DROP: Detection of RNA Outlier Pipeline")
        
        setKey(config, None, "genomeAssembly", "hg19")
        setKey(config, None, "scanBamParam", "null")
        
        # export settings
        setKey(config, None, "exportCounts", dict(), verbose=VERBOSE)
        gene_annotations = list(config["geneAnnotation"].keys())
        setKey(config, ["exportCounts"], "geneAnnotation", gene_annotations, verbose=VERBOSE)
        setKey(config, ["exportCounts"], "excludeGroups", list(), verbose=VERBOSE)
        
        # check consistency of gene annotations
        anno_incomp = set(config["exportCounts"]["geneAnnotations"]) - set(gene_annotations)
        if len(anno_incomp) > 0:
            message = f"{anno_incomp} are not valid annotation version in 'geneAnnotation'"
            message += "but required in 'exportCounts'.\n Please make sure they match."
            raise ValueError(message)
        
        if self.method is None:
            tmp_dir = submodules.getTmpDir()
        else:
            tmp_dir = submodules.getMethodPath(self.method, type_='tmp_dir')
        setKey(config, None, "tmpdir", tmp_dir)
        
        # aberrant expression
        if self.method == "AE" or self.method is None:
            setKey(config, None, "aberrantExpression", dict(), verbose=VERBOSE)
            setKey(config, ["aberrantExpression"], "fpkmCutoff", 1, verbose=VERBOSE)
            setKey(config, ["aberrantExpression"], "implementation", "autoencoder", verbose=VERBOSE)
            setKey(config, ["aberrantExpression"], "groups", None, verbose=VERBOSE)
            setKey(config, ["aberrantExpression"], "padjCutoff", .05, verbose=VERBOSE)
            setKey(config, ["aberrantExpression"], "zScoreCutoff", 0, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "maxTestedDimensionProportion", 3, verbose=VERBOSE)
        
        # aberrant splicing
        if self.method == "AS" or self.method is None:
            setKey(config, None, "aberrantSplicing", dict(), verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "groups", None, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "recount", False, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "longRead", False, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "filter", True, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "minExpressionInOneSample", 20, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "minDeltaPsi", 0, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "implementation", "PCA", verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "padjCutoff", 0.05, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "zScoreCutoff", 0.05, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "deltaPsiCutoff", 0.05, verbose=VERBOSE)
            setKey(config, ["aberrantSplicing"], "maxTestedDimensionProportion", 6, verbose=VERBOSE)
        
        # monoallelic expression
        if self.method == "MAE" or self.method is None:
            setKey(config, None, "mae", dict(), verbose=VERBOSE)
            setKey(config, ["mae"], "gatkIgnoreHeaderCheck", True, verbose=VERBOSE)
            setKey(config, ["mae"], "padjCutoff", .05, verbose=VERBOSE)
            setKey(config, ["mae"], "allelicRatioCutoff", 0.8, verbose=VERBOSE)
            setKey(config, ["mae"], "maxAF", .001, verbose=VERBOSE)
            setKey(config, ["mae"], "gnomAD", False, verbose=VERBOSE)
            setKey(config, ["mae"], "groups", None, verbose=VERBOSE)
            setKey(config, ["mae"], "qcGroups", None, verbose=VERBOSE)
        
        # commandline tools
        setKey(config, None, "tools", dict(), verbose=VERBOSE)
        setKey(config, ["tools"], "samtoolsCmd", "samtools", verbose=VERBOSE)
        setKey(config, ["tools"], "bcftoolsCmd", "bcftools", verbose=VERBOSE)
        setKey(config, ["tools"], "gatkCmd", "gatk", verbose=VERBOSE)
        
        return config
    
    def setKey(self, dict_, sub, key, default, verbose=False):
        if sub is not None:
            if not isinstance(sub, list):
                raise TypeError(f"{sub} is not of type list")
            for x in sub:
                dict_ = dict_[x]
        if key not in dict_ or dict_[key] is None:
            logger.debug(f'{key} not in config{sub}, using default')
            dict_[key] = default
        return dict_[key]


    #### FILE PARSING AND ID EXTRACTION
    
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
        
        # set ID type as string
        sample_annotation = sample_annotation.astype({
            "RNA_ID": str, "DNA_ID": str
        })
        
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
        existent = [pathlib.Path(x).exists() for x in file_mapping["FILE_PATH"]]
        if sum(existent) < file_mapping.shape[0]:
            logger.info("WARNING: there are files in the sample annotation that do not exist")
        file_mapping = file_mapping[existent].drop_duplicates()
        if file_mapping.shape[0] == 0:
            raise ValueError("No files exist in sample annotation. Please check your sample annotation.")
        
        file_mapping.to_csv(self.getProcDataDir() + "/file_mapping.csv", index=False)
        
        return file_mapping
    
    def createGroupIds(self, group_key="DROP_GROUP", file_type="RNA_BAM_FILE", sep=','):
        """
        Create a full and filtered list of RNA assay IDs subsetted by specified DROP groups
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
        grouped = {gr : df[df[group_key].str.contains(f'(^|{sep}){gr}({sep}|$)')][assay_id].tolist() 
                        for gr in groups}
        # remove groups labeled as None
        grouped = {gr : ids for gr, ids in grouped.items() if gr is not None}
        return grouped

    
    def subsetGroups(self, ids_by_group, subset_groups, warn=30, error=10):
        if subset_groups is None:
            subset = ids_by_group
        else:
            subset_groups = [subset_groups] if subset_groups.__class__ == str else subset_groups
            subset = {gr:ids for gr, ids in ids_by_group.items() if gr in subset_groups}
        
        for group in subset_groups:
            if len(subset[group]) < error:
                raise ValueError(f'Too few IDs in DROP_GROUP {group}, please ensure that it has at least {error} IDs')
            elif len(subset[group]) < warn:
                logger.info(f'WARNING: Less than {warn} IDs in DROP_GROUP {group}')
        
        return subset
    
    def createMaeIDS(self, ids_by_group, subset_groups, id_sep='--'):
        """
        create MAE IDs from smaple annotation
        """
        rna_id_by_group = self.subsetGroups(ids_by_group, subset_groups, warn=1, error=1)
        all_mae_files = self.sample_annotation[["RNA_ID", "DNA_ID"]].drop_duplicates().dropna()
        
        # subset by group
        mae_ids = {}
        for gr, rna_ids in rna_id_by_group.items():
            mae_subset = all_mae_files [all_mae_files ["RNA_ID"].isin(rna_ids)]
            vcf_rna_pairs = zip(mae_subset["DNA_ID"], mae_subset["RNA_ID"])
            mae_ids[gr] = list(map(id_sep.join, vcf_rna_pairs))
        
        return mae_ids
    

    #### GETTERS ####
    
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
            logger.debug(f"FILE NOT FOUND FOR sampleID: {sampleID} and file type {file_type}")
        return exists
    
    def getFilePath(self, sampleId, file_type):
        """
        Function for getting the file path given the sampleId and file type
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
        if self.method != 'MAE':
            self.method = 'MAE'
            self.parse()
        if not isinstance(group, str):
            group = list(group)[0]
        return self.mae_ids[group]

    def getMaeAll(self):
        """
        Get a list of all MAE IDs from the groups specified in the config.
        Useful for collecting all MAE IDs ungrouped.
        """
        if self.method != 'MAE':
            self.method = 'MAE'
            self.parse()
        all_ids = []
        groups = self.config["mae"]["groups"]
        if groups.__class__ == str:
            groups = [groups]
        for group in groups:
            all_ids.extend(self.mae_ids[group])
        return all_ids
    
    def getGeneAnnotationFile(self, annotation):
        return self.config["geneAnnotation"][annotation]
        
    def getExportGroups(self, modules=["aberrantExpression", "aberrantSplicing"]):
        groups = [] # get all active groups
        for module in modules:
          groups.extend(self.config[module]["groups"])
        exclude = self.config["exportCounts"]["excludeGroups"]
        return set(groups) - set(exclude)
        
    def getExportCountFiles(self, prefix):
        
        count_type_map = {"geneCounts":"aberrantExpression",
                          "splitCounts":"aberrantSplicing",
                          "spliceSiteOverlapCounts":"aberrantSplicing"}
        if prefix not in count_type_map.keys():
            raise ValueError(f"{prefix} not a valid file type for exported counts")
        
        datasets = self.getExportGroups([count_type_map[prefix]])
        annotations = self.config["exportCounts"]["geneAnnotations"]
        
        pattern = self.getProcResultsDir()
        if prefix == "geneCounts":
            pattern += f"/exported_counts/{{dataset}}/{prefix}_{{dataset}}--{{annotation}}.tsv.gz"
            return expand(pattern, annotation=annotations, dataset=datasets)
        else:
            pattern += f"/exported_counts/{{dataset}}/{prefix}_{{dataset}}.tsv.gz"
            return expand(pattern, dataset=datasets)

 
        

