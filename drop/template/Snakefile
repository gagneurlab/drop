from pathlib import Path
import os
import drop
import wbuild
from wbuild.createIndex import ci

projectDir = Path.cwd().resolve()
drop.checkDropVersion(projectDir)
_, projectPaths = drop.setupPaths(projectDir)
tmp_dir = projectPaths["tmpDir"]

cfg = drop.config.DropConfig(wbuild.utils.Config(), projectDir)
drop.installRPackages(cfg)
sa = cfg.sampleAnnotation
config = cfg.config_dict # legacy

conversion_dict = {
    "AberrantExpression" : "aberrant-expression-pipeline",
    "AberrantSplicing" : "aberrant-splicing-pipeline",
    "MonoallelicExpression" : "mae-pipeline",
    "rnaVariantCalling" : "rvc-pipeline"
}

def includeSnakefile(drop_module):
    drop_module.renameLocalDir()
    if drop_module.run:
        include: drop_module.getSnakefile()
        return [drop_module.name]
    else:
        return []

modules_to_run = []
modules_to_run += includeSnakefile(cfg.AE)
modules_to_run += includeSnakefile(cfg.AS)
modules_to_run += includeSnakefile(cfg.MAE)
modules_to_run += includeSnakefile(cfg.RVC)

include: drop.utils.getWBuildSnakefile()


# don't allow the dataset wildcard to contain '--' as this is a delimiter
wildcard_constraints:
    dataset="((?!--).)*"

rule all:
    input:
        expand(str(cfg.htmlOutputPath) + "/{pipeline}_index.html", pipeline=[conversion_dict[m] for m in modules_to_run]),
        rules.Index.output,
        dep_graph = cfg.get("htmlOutputPath") + "/dependency.done"

rule sampleAnnotation:
    input: rules.Pipeline_SampleAnnotation_R.output

rule exportCounts:
    input:
        cfg.exportCounts.getExportCountFiles("geneCounts"),
        cfg.exportCounts.getExportCountFiles("splicingCounts", expandPattern="k_{type}_counts", type=["j", "theta"]),
        cfg.exportCounts.getExportCountFiles("splicingCounts", expandPattern="n_{type}_counts", type=["psi5", "psi3", "theta"]),
    output:
        sampleAnnotations = cfg.exportCounts.getFiles("sampleAnnotation.tsv"),
        descriptions = cfg.exportCounts.getFiles("DESCRIPTION.txt"),
    params:
        dropConfig = cfg,
        sampleAnnotation = sa
    script: "Scripts/Pipeline/exportCountsMeta.py"

rule dependencyGraph:
    input:
        expand(str(cfg.htmlOutputPath) + "/{pipeline}_dep.svg", pipeline=[conversion_dict[m] for m in modules_to_run])
    output: touch(cfg.get("htmlOutputPath") + "/dependency.done")

rule publish_local:
    shell: "rsync -Ort {config[htmlOutputPath]} {config[webDir]}"

rule index_fa:
    input: "{path_to_ref}.fa"
    output: "{path_to_ref}.fa.fai"
    shell:
        """
        samtools faidx {input} -o {output}
        """

rule dict_fa:
    input: "{path_to_ref}.fa"
    output: "{path_to_ref}.dict"
    shell:
        """
        gatk CreateSequenceDictionary --REFERENCE {input} --OUTPUT {output}
        """

rule index_fasta:
    input: "{path_to_ref}.fasta"
    output: "{path_to_ref}.fasta.fai"
    shell:
        """
        samtools faidx {input} -o {output}
        """

rule dict_fasta:
    input: "{path_to_ref}.fasta"
    output: "{path_to_ref}.dict"
    shell:
        """
        gatk CreateSequenceDictionary --REFERENCE {input} --OUTPUT {output}
        """

