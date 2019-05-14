import pandas as pd
import os

class MyConfigParser:
    
    def __init__(self, config):
        self.config = config
        
    
    def getSampleIDs(self, experiment):
        # deprecated for all_vcf
        
        #  SAMPLE_FILE_MAPPING hast to have the following structure:
        #  [ID | FILE | ASSAY ] , ASSAY can be for example RNA_Seq
        mapping_file = self.config["SAMPLE_FILE_MAPPING"]
        df_mapping = pd.read_csv(mapping_file, sep='\t')
        
        if not list(df_mapping.columns.values)==["ID", "FILE", "ASSAY"]:
            return []
        
        # Clean and filter for experiment type e.g. EXOME_ID or RNA_seq
        df_mapping = df_mapping.dropna()
        df_mapping = df_mapping[df_mapping["ASSAY"]==experiment]
        # Check if file exists 
        df_mapping["existent"] = [os.path.exists(x) for x in df_mapping["FILE"]]
        df_mapping = df_mapping[df_mapping["existent"]]
        
        return list(df_mapping[["ID"]]) 
        
    
    def geOutriderFiles(self, experiment="RNA_ID"):
        # deprecated for outrider_files
        
        ids = self.getSampleIDs(experiment)
        
        # Get outrider Groups
        sa_file = self.config["SAMPLE_ANNOTATION"]
        df_anno = pd.read_csv(sa_file, sep='\t')
        df_anno = df_anno[df_anno["RNA_ID"].isin(ids)]
        df_outrider = df_outrider[["RNA_ID", "OUTRIDER_GROUP"]].drop_duplicates().copy()

        
        # subset by OUTRIDER_GROUP
        outrider_groups = []
        for s in set(df_outrider.OUTRIDER_GROUP):
            outrider_groups.extend(s.split(','))
        outrider_ids = {og : df_outrider.loc[anno_outrider.OUTRIDER_GROUP.str.contains('(^|,)' + og + '(,|$)'), 'RNA_ID'].tolist() for og in set(outrider_groups)}
        return outrider_ids, {og: _list for og, _list in outrider_ids.items() if len(_list) > 40}    
    
    
    def get_files_by_group(self, group):
        return expand(config["PROC_RESULTS"] + "/{{annotation}}/counts/{sampleID}.Rds", sampleID=config["outrider"][group])
    
    
    def geMaeFiles(self, rna="RNA_ID", exome="EXOME_ID"):
        # rna and exome are the names of the experiments specified in the mapping file
        
        vcf = self.getSampleIDs(exome)
        rna = self.getSampleIDs(rna)
        ### TO DO: Check if both rna and DNA for a file Exist
        return vcf, rna   

        
        
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
        
        
#def all_vcf(sa_file = config["SAMPLE_ANNOTATION"]):
#
#        anno = pd.read_csv(sa_file, sep='\t')
#
#        # subset and clean
#        anno_vcf = anno[(anno.LAB == "PROKISCH") & pd.notnull(anno.EXOME_ID)]
#        anno_vcf = anno_vcf[["EXOME_ID"]].copy()
#
#        anno_vcf['file'] = [config["RAW_DATA"] + "/" + x + "/exomicout/" for x in anno_vcf["EXOME_ID"]]
#        anno_vcf['vcf_exists'] = [os.path.exists(x) for x in anno_vcf["file"]]
#       anno_vcf = anno_vcf[anno_vcf['vcf_exists']]
#
#        return anno_vcf["EXOME_ID"].tolist()


#def outrider_files(sa_file = config["SAMPLE_ANNOTATION"]):
#        
#        anno = pd.read_csv(sa_file, sep='\t')
#        
#        # subset and clean
#        #anno_outrider = anno[(anno.LAB == "PROKISCH") & pd.notnull(anno.RNA_ID) & pd.notnull(anno.OUTRIDER_GROUP)]
#        #anno_outrider = anno_outrider[["RNA_ID", "OUTRIDER_GROUP"]].drop_duplicates().copy()
#    
#        # create filenames and ignore missing files
#        anno_outrider['file'] = [config["RAW_DATA"] + "/" + x + "/RNAout/" for x in anno_outrider["RNA_ID"]]
#        anno_outrider['file_exists'] = [os.path.exists(x) for x in anno_outrider["file"]]
#        anno_outrider = anno_outrider[anno_outrider['file_exists']]
#    
#        # subset by OUTRIDER_GROUP
#        outrider_groups = []
#       for s in set(anno_outrider.OUTRIDER_GROUP):
#            outrider_groups.extend(s.split(','))
#        outrider_ids = {og : anno_outrider.loc[anno_outrider.OUTRIDER_GROUP.str.contains('(^|,)' + og + '(,|$)'), 'RNA_ID'].tolist() for og in set(outrider_groups)}
#        return outrider_ids, {og: _list for og, _list in outrider_ids.items() if len(_list) > 40}    

        