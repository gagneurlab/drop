configfile: "make.config"
include: ".wBuild/wBuild.snakefile"

rule all:
	input: rules.Index.output, "Output/html/readme.html"
