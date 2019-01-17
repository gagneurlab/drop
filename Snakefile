configfile: "wbuild.yaml"
include: ".wBuild/wBuild.snakefile"

htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")

rule count:
    input: expand(config["PROC_RESULTS"] + "/{annotation}/counts/total_counts_{strand}.Rds", annotation=config["ANNOTATIONS"], strand=['ss', 'ns'])
