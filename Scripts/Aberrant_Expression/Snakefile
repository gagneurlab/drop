configfile: "wbuild.yaml"
include: ".wBuild/wBuild.snakefile"

htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

rule all:
	input: rules.Index.output, htmlOutputPath + "/readme.html"
	output: touch("Output/all.done")
