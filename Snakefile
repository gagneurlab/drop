import pandas as pd
import os
import numpy as np

configfile: "wbuild.yaml"

# sample annotation
subworkflow standardFileNames:
    workdir:
        "../sample_annotation"
    snakefile:
        "../sample_annotation/Snakefile"
    configfile:
        "../sample_annotation/wbuild.yaml"

# 1 aberrant Expression
subworkflow aberrantExp:
    workdir:
        "../aberrant-expression-pipeline"
    snakefile:
        "../aberrant-expression-pipeline/Snakefile"
    configfile:
        "../aberrant-expression-pipeline/wbuild.yaml"

# 2 aberrant Splicing
subworkflow aberrantSplicing:
    workdir:
        "../aberrant-splicing-pipeline"
    snakefile:
        "../aberrant-splicing-pipeline/Snakefile"
    configfile:
        "../aberrant-splicing-pipeline/wbuild.yaml"

# 3 mae
subworkflow mae:
    workdir:
        "../mae-pipeline"
    snakefile:
        "../mae-pipeline/Snakefile"
    configfile:
        "../mae-pipeline/wbuild.yaml"

# 4 variants
subworkflow variants:
    workdir:
        "../variant-annotation-pipeline"
    snakefile:
        "../variant-annotation-pipeline/Snakefile"
    configfile:
        "../variant-annotation-pipeline/wbuild.yaml"

# 5 proteomics
subworkflow proteomics:
    workdir:
        "../proteomics-pipeline"
    snakefile:
        "../proteomics-pipeline/Snakefile"
    configfile:
        "../proteomics-pipeline/wbuild.yaml"





include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
htmlOutputPath = "Output/html" # config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"


rule all:
    input: 
        rules.Index.output, # rule.Index.output is  "Output/html/index.html"
        htmlOutputPath + "/readme.html",
        variants( "../variant-annotation-pipeline/Output/html/index.html"),
        variants("../variant-annotation-pipeline/Output/html/readme.html"),
        aberrantExp( "../aberrant-expression-pipeline/Output/html/index.html"),
        aberrantExp("../aberrant-expression-pipeline/Output/html/readme.html"),
        aberrantSplicing( "../aberrant-splicing-pipeline/Output/html/index.html"),
        aberrantSplicing("../aberrant-splicing-pipeline/Output/html/readme.html"),
        mae( "../mae-pipeline/Output/html/index.html"),
        mae("../mae-pipeline/Output/html/readme.html"),
        proteomics( "../proteomics-pipeline/Output/html/index.html"),
        proteomics("../proteomics-pipeline/Output/html/readme.html")
    output: 
        touch("Output/all.done"),
        touch("../variant-annotation-pipeline/Output/all.done"),
        touch("../aberrant-expression-pipeline/Output/all.done"),
        touch("../aberrant-splicing-pipeline/Output/all.done"),
        touch("../mae-pipeline/Output/all.done"),
        touch("../proteomics-pipeline/Output/all.done")
        

rule variants:
    input:
        variants( "../variant-annotation-pipeline/Output/html/index.html"),
        variants("../variant-annotation-pipeline/Output/html/readme.html")
    output:
        touch("../variant-annotation-pipeline/Output/all.done")
     
        
rule aberrant_expression:
    input:
        aberrantExp("../aberrant-expression-pipeline/Output/html/index.html"),
        aberrantExp("../aberrant-expression-pipeline/Output/html/readme.html")
    output:
        touch("../aberrant-expression-pipeline/Output/all.done")
        
rule aberrant_splicing:
    input:
        aberrantSplicing( "../aberrant-splicing-pipeline/Output/html/index.html"),
        aberrantSplicing("../aberrant-splicing-pipeline/Output/html/readme.html")
    output:
        touch("../aberrant-splicing-pipeline/Output/all.done")

rule mae:
    input: 
        mae( "../mae-pipeline/Output/html/index.html"),
        mae("../mae-pipeline/Output/html/readme.html")
    output:
        touch("../mae-pipeline/Output/all.done")
        
rule proteomics:
    input:
        proteomics( "../proteomics-pipeline/Output/html/index.html"),
        proteomics("../proteomics-pipeline/Output/html/readme.html")
    output:
        touch("../proteomics-pipeline/Output/all.done")
        
