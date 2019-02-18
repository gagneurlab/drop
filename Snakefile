configfile: "wbuild.yaml"
include: ".wBuild/wBuild.snakefile"

import pandas as pd
import os
import numpy as np

htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"


rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html", "Output/mae.done"
    output: touch("Output/all.done")

rule count:
    input: expand(config["PROC_RESULTS"] + "/{annotation}/counts/total_counts_{strand}.Rds", annotation=config["ANNOTATIONS"], strand=['ss', 'ns'])

def mae_files(wildcards):
    anno = pd.read_table(config["SAMPLE_ANNOTATION"])

    # subset and clean
    anno_mae = anno[anno["LAB"] == "PROKISCH"]
    anno_mae = anno_mae[pd.notnull(anno_mae.EXOME_ID)]
    anno_mae = anno_mae[pd.notnull(anno_mae.RNA_ID)]
    anno_mae = anno_mae[["EXOME_ID", "RNA_ID"]]

    # create file names
    anno_mae['rna_file'] = [config["RAW_DATA"] + "/" + x + "/RNAout/paired-endout/stdFilenames/" + x + ".bam" for x in anno_mae["RNA_ID"]]
    anno_mae['vcf_file'] = [config["RAW_DATA"] + "/" + x + "/exomicout/paired-endout/stdFilenames/" + x + ".vcf.gz" for x in anno_mae["EXOME_ID"]]

    # check for missing files
    anno_mae['vcf_exists'] = [os.path.exists(x) for x in anno_mae["vcf_file"]]
    anno_mae['rna_exists'] = [os.path.exists(x) for x in anno_mae["rna_file"]]
    anno_mae = anno_mae[anno_mae['vcf_exists'] & anno_mae['rna_exists']]
    
    out_files = config["PROC_DATA"] + '/mae/' + anno_mae["EXOME_ID"] + '-' + anno_mae["RNA_ID"] + '.Rds'
    
    return out_files.tolist()
    
rule mae:
    input: mae_files
    output: touch("Output/mae.done")

