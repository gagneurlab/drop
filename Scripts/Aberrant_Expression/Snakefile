import pandas as pd
import os
import numpy as np

configfile: "wbuild.yaml"

subworkflow standardFileNames:
    workdir:
        "../../../sample_annotation"
    snakefile:
        "../../../sample_annotation/Snakefile"
    configfile:
        "../../../sample_annotation/wbuild.yaml"


def outrider_files(sa_file = config["SAMPLE_ANNOTATION"]):

    anno = pd.read_csv(sa_file, sep='\t')
    
    # subset and clean
    anno_outrider = anno[(anno.LAB == "PROKISCH") & pd.notnull(anno.RNA_ID) & pd.notnull(anno.OUTRIDER_GROUP)]
    anno_outrider = anno_outrider[["RNA_ID", "OUTRIDER_GROUP"]].drop_duplicates().copy()

    # create filenames and ignore missing files
    anno_outrider['file'] = [config["RAW_DATA"] + "/" + x + "/RNAout/" for x in anno_outrider["RNA_ID"]]
    anno_outrider['file_exists'] = [os.path.exists(x) for x in anno_outrider["file"]]
    anno_outrider = anno_outrider[anno_outrider['file_exists']]

    # subset by OUTRIDER_GROUP
    outrider_groups = []
    for s in set(anno_outrider.OUTRIDER_GROUP):
        outrider_groups.extend(s.split(','))
    outrider_ids = {og : anno_outrider.loc[anno_outrider.OUTRIDER_GROUP.str.contains('(^|,)' + og + '(,|$)'), 'RNA_ID'].tolist() for og in set(outrider_groups)}
    return outrider_ids, {og: _list for og, _list in outrider_ids.items() if len(_list) > 40}    

def get_files_by_group(group):
    return expand(config["PROC_RESULTS"] + "/{{annotation}}/counts/{sampleID}.Rds", sampleID=config["outrider"][group])


def all_vcf(sa_file = config["SAMPLE_ANNOTATION"]):

    anno = pd.read_csv(sa_file, sep='\t')
    
    # subset and clean
    anno_vcf = anno[(anno.LAB == "PROKISCH") & pd.notnull(anno.EXOME_ID)]
    anno_vcf = anno_vcf[["EXOME_ID"]].copy()

    anno_vcf['file'] = [config["RAW_DATA"] + "/" + x + "/exomicout/" for x in anno_vcf["EXOME_ID"]]
    anno_vcf['vcf_exists'] = [os.path.exists(x) for x in anno_vcf["file"]]
    anno_vcf = anno_vcf[anno_vcf['vcf_exists']]
    
    return anno_vcf["EXOME_ID"].tolist()
    
    
def mae_files(sa_file = config["SAMPLE_ANNOTATION"]):
    
    anno = pd.read_csv(sa_file, sep='\t')
    
    # subset and clean
    anno_mae = anno[anno["LAB"] == "PROKISCH"]
    anno_mae = anno_mae[pd.notnull(anno_mae.EXOME_ID)]
    anno_mae = anno_mae[pd.notnull(anno_mae.RNA_ID)]
    anno_mae = anno_mae[["EXOME_ID", "RNA_ID"]].copy()

    # create file names
    # anno_mae['rna_file'] = [config["RAW_DATA"] + "/" + x + "/RNAout/paired-endout/stdFilenames/" + x + ".bam" for x in anno_mae["RNA_ID"]]
    # anno_mae['vcf_file'] = [config["RAW_DATA"] + "/" + x + "/exomicout/paired-endout/stdFilenames/" + x + ".vcf.gz" for x in anno_mae["EXOME_ID"]]
    
    anno_mae['rna_file'] = [config["RAW_DATA"] + "/" + x + "/RNAout" for x in anno_mae["RNA_ID"]]
    anno_mae['vcf_file'] = [config["RAW_DATA"] + "/" + x + "/exomicout" for x in anno_mae["EXOME_ID"]]

    # check for missing files
    anno_mae['vcf_exists'] = [os.path.exists(x) for x in anno_mae["vcf_file"]]
    anno_mae['rna_exists'] = [os.path.exists(x) for x in anno_mae["rna_file"]]
    anno_mae = anno_mae[anno_mae['vcf_exists'] & anno_mae['rna_exists']]
    
    vcf = anno_mae["EXOME_ID"] 
    rna = anno_mae["RNA_ID"]
    
    return vcf.tolist(), rna.tolist()


# set config variables
#mae
vcfs, rnas = mae_files()
config["vcfs"] = vcfs
config["rnas"] = rnas
config["mae_ids"] = list(map('-'.join, zip(vcfs, rnas)))

#outrider
outrider_all_ids, outrider_filtered = outrider_files()
config["outrider"] = outrider_all_ids
config["outrider_filtered"] = outrider_filtered

include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
#htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
htmlOutputPath = "Output/html"


rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")
