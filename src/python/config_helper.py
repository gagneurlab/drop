import pandas as pd
import os
import wbuild.utils as wbu

class ConfigHelper:
    
    def __init__(self, config):

        if not config:
            wconf = wbu.Config()
            config = wconf.conf_dict
        self.config = config
        
        # sample-file mappping: reading and cleaning 
        #  SAMPLE_FILE_MAPPING has to have the following structure:
        #  [ID | FILE | ASSAY ] , ASSAY can be for example RNA_Seq
        df_mapping = pd.read_csv(self.config["SAMPLE_FILE_MAPPING"], sep='\t')
        if not list(df_mapping.columns.values)==["ID", "FILE", "ASSAY"]:
            print("File does not correspond to required format with columns [ID | FILE | ASSAY]")
        
        df_mapping = df_mapping.dropna()
        
        # Check if file exists 
        df_mapping["existent"] = [os.path.exists(x) for x in df_mapping["FILE"]]
        df_mapping = df_mapping[df_mapping["existent"]]
        self.sample_file_mapping = df_mapping
        #print(self.sample_file_mapping.tail())
        
        # sample annotation
        #  SAMPLE_ANNOTATION must have assay names as specified in sample-file mappping for ID columns
        sa_file = self.config["SAMPLE_ANNOTATION"]
        self.sample_annotation = pd.read_csv(sa_file, sep='\t')
        
        # OUTRIDER ids
        self.outrider_all, self.outrider_filtered = self.createOutriderIds(min_ids=self.config["min_outrider_ids"])
    
        #print(self.sample_file_mapping.head(5))
    
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
        # deprecated for all_vcf
        return list(self.sample_file_mapping[self.sample_file_mapping["ASSAY"] == experiment]["ID"]) 
    
    
    def checkFileExists(self, sampleId, assay):
        #print(sampleId, assay)
        x = self.sample_file_mapping[(self.sample_file_mapping["ASSAY"]==assay) & ((self.sample_file_mapping["ID"]==sampleId))]["FILE"]#.iloc[0]
        if len(x)<1:
          print("ENTRY NOT FOUND in sample_file_mapping for sampleId: {} and assay: {}".format(sampleId, assay))
          return False
        exists = os.path.exists(x.iloc[0])
        if not exists:
          print("FILE NOT FOUND FOR sampleID: ", sampleId, "and assay", assay)
        return exists
        
    """
    Returns vcf and rna files for MAE pipeline
    """
    def getMaeIDs(self):
        # rna and exome are the names of the experiments specified in the mapping file
        rna_assay = self.config["rna_assay"]
        dna_assay = self.config["dna_assay"]
        #print("rna_assay: <{}>, dna_assay: <{}>".format(rna_assay,dna_assay))
        # return nothing, if there aren't any exomes
        
        self.sample_annotation = self.sample_annotation[pd.notnull(self.sample_annotation[rna_assay])]
        self.sample_annotation = self.sample_annotation[pd.notnull(self.sample_annotation[dna_assay])]
        
        self.sample_annotation['vcf_exists'] = [self.checkFileExists(x, dna_assay) for x in list(self.sample_annotation[dna_assay])]
        self.sample_annotation['rna_exists'] = [self.checkFileExists(x, rna_assay) for x in list(self.sample_annotation[rna_assay])]
        self.sample_annotation = self.sample_annotation[self.sample_annotation['vcf_exists'] & self.sample_annotation['rna_exists']]

        vcfs = list(self.sample_annotation[dna_assay]) 
        rnas = list(self.sample_annotation[rna_assay])
        print("length of vcfs: {}, and rnas: {}".format(len(vcfs), len(rnas)))
        return vcfs, rnas 

      
        
    
    """Function for getting the file path given the sampleId and assay
    @param sampleId: ID of sample
    @param assay: either "rna_assay", "dna_assay", as specified in the config
    """
    def getFilePath(self, sampleId, assay_name):
        #print("In function getFilePath: --sampleID: {}, assay_name: {}".format(sampleId, assay_name))
        assay = self.config[assay_name]
        #deprecated for stdFileNames from subworkflow sample_annotation
        self.sample_file_mapping["ID"] = self.sample_file_mapping["ID"].astype(str)
        
        #print(self.sample_file_mapping)
        path = (self.sample_file_mapping[(self.sample_file_mapping["ASSAY"] == assay) & (self.sample_file_mapping["ID"] == sampleId)]["FILE"]).iloc[0]
        
        #print("file path", path)
        return path
      
    
    """
    Create a full and filtered list of RNA assay IDs subsetted by specified OUTRIDER groups
    """
    def createOutriderIds(self, min_ids=40):
        # deprecated for outrider_files
        
        outrider_group_col = self.config["outrider_group"]
        rna_assay = self.config["rna_assay"]
        ids = self.getSampleIDs(self.config["rna_assay"])
        
        # Get unique outrider Groups
        df_outrider = self.sample_annotation[self.sample_annotation[rna_assay].isin(ids)]
        df_outrider = df_outrider[[rna_assay, outrider_group_col]].drop_duplicates().copy()
        #print("df_outrider", df_outrider.head(5))
        
        # assumes that OUTRIDER groups are comma-separated
        # get unique group names
        outrider_groups = []
        for s in set(df_outrider[outrider_group_col]):
            outrider_groups.extend(s.split(','))
        #print("outrider_groups", outrider_groups)
            
        # collect IDs per group
        outrider_ids = {og : df_outrider.loc[
            df_outrider[outrider_group_col].str.contains('(^|,)' + og + '(,|$)'),
            rna_assay].tolist() for og in set(outrider_groups)}
        outrider_filtered = {og: _list for og, _list in outrider_ids.items() if len(_list) > min_ids}
        
        #print("outrider ids", outrider_ids)
        #print()
        #print("outrider filtered", outrider_filtered)
        return outrider_ids, outrider_filtered
    
    """
    Get lists of IDs per OUTRIDER group as (all ids, filtered ids)
    """
    def getOutriderIds(self):
        return self.outrider_all, self.outrider_filtered
    
    """ Wrapper for getting all count files for the specified OUTRIDER group
    
    @param annotation: wildcard specifying the gene annotation version from config["GENE_ANNOTATION_NAMES"]
    @param group: group label of sample, under outrider column in sample annotation
    """
    def getCountFileByOutriderGroup(self, annotation, group):
        res = []
        for sampleid in self.outrider_all[group]:
            res.append(f"{self.getProcDataDir()}/{annotation}/counts/{str(sampleid)}.Rds")
        return res
    
        #return expand(self.getProcResultsDir() + "/{{annotation}}/counts/{sampleID}.Rds", sampleID=self.outrider_all[group])
        # DOES NOT work because expand is defined in Snakemake and not in Python
    
    def getGeneAnnotationFile(self, annotation):
        i = self.config["GENE_ANNOTATION_NAMES"].index(annotation)
        return self.config["GENE_ANNOTATION"][i]
        
    def getGeneInfoFile(self, annotation):
        i = self.config["GENE_ANNOTATION_NAMES"].index(annotation)
        return self.config["GENE_INFO"][i]
        
    def getCountRangesFile(self, annotation):
        i = self.config["GENE_ANNOTATION_NAMES"].index(annotation)
        return self.config["COUNT_RANGES"][i]
    

        
        
#    def mae_files(sa_file = config["SAMPLE_ANNOTATION"]):
#        
#        anno = pd.read_csv(sa_file, sep='\t')
#        
#        # subset and clean
#        anno_mae = anno[anno["LAB"] == "PROKISCH"]
#        anno_mae = anno_mae[pd.notnull(anno_mae.EXOME_ID)]
#        anno_mae = anno_mae[pd.notnull(anno_mae.RNA_ID)]
#        anno_mae = anno_mae[["EXOME_ID", "RNA_ID"]].copy()
#    
#        # create file names
#        # anno_mae['rna_file'] = [config["RAW_DATA"] + "/" + x + "/RNAout/paired-endout/stdFilenames/" + x + ".bam" for x in anno_mae["RNA_ID"]]
#        # anno_mae['vcf_file'] = [config["RAW_DATA"] + "/" + x + "/exomicout/paired-endout/stdFilenames/" + x + ".vcf.gz" for x in anno_mae["EXOME_ID"]]
#        
#       anno_mae['rna_file'] = [config["RAW_DATA"] + "/" + x + "/RNAout" for x in anno_mae["RNA_ID"]]
#       anno_mae['vcf_file'] = [config["RAW_DATA"] + "/" + x + "/exomicout" for x in anno_mae["EXOME_ID"]]
#    
#        # check for missing files
#        anno_mae['vcf_exists'] = [os.path.exists(x) for x in anno_mae["vcf_file"]]
#        anno_mae['rna_exists'] = [os.path.exists(x) for x in anno_mae["rna_file"]]
#        anno_mae = anno_mae[anno_mae['vcf_exists'] & anno_mae['rna_exists']]
#        
#        vcf = anno_mae["EXOME_ID"] 
#        rna = anno_mae["RNA_ID"]
#        
#        return vcf.tolist(), rna.tolist()


# TO DO: Check if BOTH rna and dna files exist
#for i in range(len(rnas)):
#    for i in range(len(vcfs)):
#        
#        vcf_exists = os.path.exists(self.sample_file_mapping[(self.sample_file_mapping["ID"]==vcfs[i] & (self.sample_file_mapping["ASSAY"]==dna_assay)]["FILE"])
#        rna_exists = os.path.exists(self.sample_file_mapping[(self.sample_file_mapping["ID"]==rnas[i] & (self.sample_file_mapping["ASSAY"]==rna_assay)]["FILE"])
#        
#        if not vcf_exists:
#            print("Missing vcf File for sampleID", vcfs[i])
#        if not rna_exists:
#            print("Missing rna File for sampleID", rnas[i])
#            
#        if not (vcf_exists and rna_exists):
#            rnas.pop(i)
#            vcfs.pop(i)
    
 
        
