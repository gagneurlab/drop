### SNAKEFILE MONOALLELIC EXPRESSION
import os
import drop
import pathlib

METHOD = 'MAE'
SCRIPT_ROOT = drop.getMethodPath(METHOD, type_='workdir', str_=False)
CONF_FILE = drop.getConfFile()

parser = drop.config(config, METHOD)
config = parser.parse()
include: config['wBuildPath'] + "/wBuild.snakefile"

###### FUNCTIONS ######
def fasta_dict(fasta_file):
    return fasta_file.split('.')[0] + ".dict"

def getVcf(rna_id, vcf_id="qc"):
    if vcf_id == "qc":
        return config["mae"]["qcVcf"]
    else:
        return parser.getProcDataDir() + f"/mae/snvs/{vcf_id}--{rna_id}.vcf.gz"
        
def getQC(format):
    if format == "UCSC":
        return config["mae"]["qcVcf"]
    elif format == "NCBI":
        return parser.getProcDataDir() + "/mae/qc_vcf_ncbi.vcf.gz"
    else:
        raise ValueError(f"getQC: {format} is an invalid chromosome format")

def getChrMap(SCRIPT_ROOT, conversion):
    if conversion == 'ncbi2ucsc':
        return SCRIPT_ROOT/"resource"/"chr_NCBI_UCSC.txt"
    elif conversion == 'ucsc2ncbi':
        return SCRIPT_ROOT/"resource"/"chr_UCSC_NCBI.txt"
    else:
        raise ValueError(f"getChrMap: {conversion} is an invalid conversion option")
        
def getScript(type, name):
    return SCRIPT_ROOT/"Scripts"/type/name
######

rule all:
    input: 
        rules.Index.output, 
        config["htmlOutputPath"] + "/mae_readme.html",
        rules.Scripts_MAE_Datasets_R.output,
        rules.Scripts_QC_Datasets_R.output
    output: touch(drop.getMethodPath(METHOD, type_='final_file'))

rule sampleQC:
    input: rules.Scripts_QC_Datasets_R.output
    output: touch(drop.getTmpDir() + "/sampleQC.done")

rule create_dict:
    input: config['mae']['genome']
    output: fasta_dict(config['mae']['genome'])
    shell: "gatk CreateSequenceDictionary --REFERENCE {input[0]}"
        
## MAE
rule create_SNVs:
    input:
        ncbi2ucsc = getChrMap(SCRIPT_ROOT, "ncbi2ucsc"),
        ucsc2ncbi = getChrMap(SCRIPT_ROOT, "ucsc2ncbi"),
        vcf_file  = lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, 
                    file_type='DNA_VCF_FILE'),
        bam_file  = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, 
                    file_type='RNA_BAM_FILE'),
        script    = getScript("MAE", "filterSNVs.sh")
    output:
        snvs_filename=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz",
        snvs_index=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz.tbi"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} {input.vcf_file} \
        {wildcards.vcf} {input.bam_file} {output.snvs_filename} \
        {config[tools][bcftoolsCmd]} {config[tools][samtoolsCmd]}
        """

rule allelic_counts: 
    input:
        ncbi2ucsc = getChrMap(SCRIPT_ROOT, "ncbi2ucsc"),
        ucsc2ncbi = getChrMap(SCRIPT_ROOT, "ucsc2ncbi"),
        vcf_file  = lambda wildcards: getVcf(wildcards.rna, wildcards.vcf),
        bam_file  = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, 
                    file_type='RNA_BAM_FILE'),
        fasta     = config['mae']['genome'],
        dict      = fasta_dict(config['mae']['genome']),
        script    = getScript("MAE", "ASEReadCounter.sh")
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file} {input.bam_file} {wildcards.vcf}--{wildcards.rna} \
        {input.fasta} {config[mae][gatkIgnoreHeaderCheck]} {output.counted} \
        {config[tools][bcftoolsCmd]}
        """
## QC
rule renameChrQC:
    input:
        ucsc2ncbi = getChrMap(SCRIPT_ROOT, "ucsc2ncbi"),
        ncbi_vcf = getQC(format="UCSC")
    output:
        ncbi_vcf = getQC(format="NCBI")
    shell:
        """
        bcftools={config[tools][bcftoolsCmd]}
        echo 'converting from UCSC to NCBI format'
        $bcftools annotate --rename-chrs {input.ucsc2ncbi} {input.ncbi_vcf} \
            | bgzip > {output.ncbi_vcf}
        $bcftools index -t {output.ncbi_vcf}
        """

rule allelic_counts_qc: 
    input:
        ncbi2ucsc = getChrMap(SCRIPT_ROOT, "ncbi2ucsc"),
        ucsc2ncbi = getChrMap(SCRIPT_ROOT, "ucsc2ncbi"),
        vcf_file_ucsc = getQC(format="UCSC"),
        vcf_file_ncbi = getQC(format="NCBI"),
        bam_file      = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, 
                        file_type='RNA_BAM_FILE'),
        fasta         = config['mae']['genome'],
        dict          = fasta_dict(config['mae']['genome']),
        script_qc = getScript("QC", "ASEReadCounter.sh"),
        script_mae = getScript("MAE", "ASEReadCounter.sh")
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts/qc_{rna}.csv.gz"
    shell:
        """
        {input.script_qc} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file_ucsc} {input.vcf_file_ncbi} {input.bam_file} \
        {wildcards.rna} {input.fasta} {config[mae][gatkIgnoreHeaderCheck]} \
        {output.counted} {config[tools][bcftoolsCmd]} \
        {config[tools][samtoolsCmd]} {input.script_mae}
        """

####
rulegraph_filename = f'{config["htmlOutputPath"]}/{METHOD}_rulegraph'
rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        svg = f"{rulegraph_filename}.svg",
        png = f"{rulegraph_filename}.png"
    shell:
        """
        snakemake --configfile {CONF_FILE} --rulegraph | dot -Tsvg > {output.svg}
        snakemake --configfile {CONF_FILE} --rulegraph | dot -Tpng > {output.png}
        """

rule unlock:
    output: touch(drop.getMethodPath(METHOD, type_="unlock"))
    shell: "snakemake --unlock --configfile {CONF_FILE}"

