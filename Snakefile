### SNAKEFILE GENETIC DIAGNOSIS PIPELINE
import os
from config_parser import ConfigHelper

## ADD tmp/ DIR
if not os.path.exists('tmp'):
    os.makedirs('tmp')

### Write one config file for every subworkflow with a diferent index name
import oyaml
with open('tmp/config_aberrant_expression.yaml', 'w') as yaml_file:
    config_ae = config.copy()
    oyaml.dump(config_ae, yaml_file, default_flow_style=False)
    
with open('tmp/config_aberrant_splicing.yaml', 'w') as yaml_file:
    config_as = config.copy()
    oyaml.dump(config_as, yaml_file, default_flow_style=False)
    
with open('tmp/config_mae.yaml', 'w') as yaml_file:
    config_mae = config.copy()
    oyaml.dump(config_mae, yaml_file, default_flow_style=False)    


print("In Snakefile from genetic_diagnosis",config)
parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile

htmlOutputPath = config["htmlOutputPath"]

# 1 aberrant Expression
subworkflow aberrantExp:
    workdir:
        "submodules/aberrant-expression-pipeline"
    snakefile:
        "submodules/aberrant-expression-pipeline/Snakefile"
    configfile:
        "tmp/config_aberrant_expression.yaml"

# 2 aberrant Splicing
subworkflow aberrantSplicing:
    workdir:
        "submodules/aberrant-splicing-pipeline"
    snakefile:
        "submodules/aberrant-splicing-pipeline/Snakefile"
    configfile:
        "tmp/config_aberrant_splicing.yaml"

# 3 mae
subworkflow mae:
    workdir:
        "submodules/mae-pipeline"
    snakefile:
        "submodules/mae-pipeline/Snakefile"
    configfile:
        "tmp/config_mae.yaml"

include: os.getcwd() + "/.wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
include: os.getcwd() + "/.wBuild/wBuild.snakefile" 
if not os.path.exists('tmp'):
    os.makedirs('tmp')

rule all:
    input: 
        aberrantExp("tmp/aberrant_expression.done"),
        aberrantSplicing("tmp/aberrant_splicing.done"),
        mae("tmp/mae.done")        
    output:
        touch("tmp/gdp_all.done")

######## Do not delete this!!!
rule aberrant_expression:
    input:
        aberrantExp("tmp/aberrant_expression.done")
    output:
        touch("tmp/gdp_aberrant_expression.done")
        
rule aberrant_splicing:
    input:
        aberrantSplicing("tmp/aberrant_splicing.done")
    output:
        touch("tmp/gdp_aberrant_splicing.done")

rule MAE:
    input:
        mae("tmp/mae.done")
    output:
        touch("tmp/gdp_mae.done")
        

rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

