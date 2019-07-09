import os
from config_parser import ConfigHelper

parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
htmlOutputPath = config["htmlOutputPath"]

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


include: os.getcwd() + "/.wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
if not os.path.exists('tmp'):
    os.makedirs('tmp')

rule all:
    input: 
        aberrantExp(htmlOutputPath  + "/aberrant_expression.done"),
        aberrantSplicing(htmlOutputPath  + "/aberrant_splicing.done"),
        mae(htmlOutputPath  + "/mae.done")
    output:
        touch(htmlOutputPath + "/../all.done"),



rule aberrant_expression:
    input:
        aberrantExp(htmlOutputPath + "/readme.html")
    output:
        touch(htmlOutputPath  + "/aberrant_expression.done")
        
rule aberrant_splicing:
    input:
        aberrantSplicing(htmlOutputPath + "/readme.html")
    output:
        touch(htmlOutputPath + "/aberrant_splicing.done")

rule mae:
    input: 
        mae(htmlOutputPath + "/readme.html")
    output:
        touch(htmlOutputPath + "/mae.done")
        

rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

