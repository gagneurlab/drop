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
htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"


rule all:
    input: 
        rules.Index.output, # rule.Index.output is  "Output/html/index.html"
        htmlOutputPath + "/readme.html",
        aberrantExp(htmlOutputPath + "/index.html"),
        aberrantExp(htmlOutputPath + "/readme.html"),
        aberrantSplicing(htmlOutputPath + "/index.html"),
        aberrantSplicing(htmlOutputPath + "/readme.html"),
        mae(htmlOutputPath + "/index.html"),
        mae(htmlOutputPath + "/readme.html")
    output: 
        touch(htmlOutputPath + "/../all.done"),
        touch("submodules/aberrant-expression-pipeline/" + htmlOutputPath + "/../all.done"),
        touch("submodules/aberrant-splicing-pipeline/" + htmlOutputPath + "/../all.done"),
        touch("submodules/mae-pipeline/" + htmlOutputPath + "/../all.done")



rule aberrant_expression:
    input:
        aberrantExp(htmlOutputPath + "/index.html"),
        aberrantExp(htmlOutputPath + "/readme.html")
    output:
        touch("submodules/aberrant-expression-pipeline/" + htmlOutputPath  + "/../all.done")
        
rule aberrant_splicing:
    input:
        aberrantSplicing(htmlOutputPath + "/index.html"),
        aberrantSplicing(htmlOutputPath + "/readme.html")
    output:
        touch("submodules/aberrant-splicing-pipeline/" + htmlOutputPath + "/../all.done")

rule mae:
    input: 
        mae(htmlOutputPath + "/index.html"),
        mae(htmlOutputPath + "/readme.html")
    output:
        touch("submodules/mae-pipeline/" + htmlOutputPath + "/../all.done")
        

rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

