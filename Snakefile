import pandas as pd
import os
import numpy as np


# 1 aberrant Expression
subworkflow aberrantExp:
    workdir:
        "submodules/aberrant-expression-pipeline"
    snakefile:
        "submodules/aberrant-expression-pipeline/Snakefile"

# 2 aberrant Splicing
subworkflow aberrantSplicing:
    workdir:
        "submodules/aberrant-splicing-pipeline"
    snakefile:
        "submodules/aberrant-splicing-pipeline/Snakefile"

# 3 mae
subworkflow mae:
    workdir:
        "submodules/mae-pipeline"
    snakefile:
        "submodules/mae-pipeline/Snakefile"

## 4 variants
#subworkflow variants:
#    workdir:
#        "submodules/variant-annotation-pipeline"
#    snakefile:
#        "submodules/variant-annotation-pipeline/Snakefile"



include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
htmlOutputPath = "Output/html" # config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"


rule all:
    input: 
        rules.Index.output, # rule.Index.output is  "Output/html/index.html"
        htmlOutputPath + "/readme.html",
        aberrantExp( "submodules/aberrant-expression-pipeline/Output/html/index.html"),
        aberrantExp("submodules/aberrant-expression-pipeline/Output/html/readme.html"),
        aberrantSplicing( "submodules/aberrant-splicing-pipeline/Output/html/index.html"),
        aberrantSplicing("submodules/aberrant-splicing-pipeline/Output/html/readme.html"),
        mae( "submodules/mae-pipeline/Output/html/index.html"),
        mae("submodules/mae-pipeline/Output/html/readme.html")
    output: 
        touch("Output/all.done"),
        touch("submodules/aberrant-expression-pipeline/Output/all.done"),
        touch("submodules/aberrant-splicing-pipeline/Output/all.done"),
        touch("submodules/mae-pipeline/Output/all.done")



rule aberrant_expression:
    input:
        aberrantExp("submodules/aberrant-expression-pipeline/Output/html/index.html"),
        aberrantExp("submodules/aberrant-expression-pipeline/Output/html/readme.html")
    output:
        touch("submodules/aberrant-expression-pipeline/Output/all.done")
        
rule aberrant_splicing:
    input:
        aberrantSplicing( "submodules/aberrant-splicing-pipeline/Output/html/index.html"),
        aberrantSplicing("submodules/aberrant-splicing-pipeline/Output/html/readme.html")
    output:
        touch("submodules/aberrant-splicing-pipeline/Output/all.done")

rule mae:
    input: 
        mae( "submodules/mae-pipeline/Output/html/index.html"),
        mae("submodules/mae-pipeline/Output/html/readme.html")
    output:
        touch("submodules/mae-pipeline/Output/all.done")
        
rule proteomics:
    input:
        proteomics( "submodules/proteomics-pipeline/Output/html/index.html"),
        proteomics("submodules/proteomics-pipeline/Output/html/readme.html")
    output:
        touch("submodules/proteomics-pipeline/Output/all.done")

rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

