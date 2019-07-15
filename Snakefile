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


include: os.getcwd() + "/.wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
if not os.path.exists('tmp'):
    os.makedirs('tmp')

rule all:
    input: 
        aberrantExp(rules.Index.output),
        aberrantSplicing(rules.Index.output),
        mae(rules.Index.output)        
    output:
        touch(htmlOutputPath + "/../all.done"),


######## Do not delete this!!!
rule aberrant_expression:
    input:
        aberrantExp(rules.Index.output)
    output:
        touch(htmlOutputPath  + "/aberrant_expression.done")
        
rule aberrant_splicing:
    input:
        aberrantSplicing(rules.Index.output)
    output:
        touch(htmlOutputPath + "/aberrant_splicing.done")

rule mae:
    input: 
        mae(rules.Index.output)
    output:
        touch(htmlOutputPath + "/mae.done")
        

rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

