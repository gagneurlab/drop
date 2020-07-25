WORKDIR = cfg.AE.getWorkdir(str_=False)

rule aberrantExpression:
    input:
        expand(
            cfg.getHtmlFromScript(WORKDIR / "Counting" / "Datasets.R"),
            annotation=cfg.getGeneVersions()
        ),
        expand(
            cfg.getProcessedResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv",
                annotation=cfg.getGeneVersions(), dataset=cfg.AE.groups
        ),
        expand(
            cfg.getHtmlFromScript(WORKDIR / "OUTRIDER" / "Datasets.R"),
            annotation=cfg.getGeneVersions()
        )

rule aberrantExpression_bamStats:
    input: 
        bam = lambda wildcards: sa.getFilePath(wildcards.sampleID, "RNA_BAM_FILE"),
        ucsc2ncbi = WORKDIR / "resource" / "chr_UCSC_NCBI.txt"
    output:
        cfg.processedDataDir / "aberrant_expression" / "{annotation}" / "coverage" / "{sampleID}.tsv"
    params:
        samtools = config["tools"]["samtoolsCmd"]
    shell:
        """
        # identify chromosome format
        if {params.samtools} idxstats {input.bam} | grep -qP "^chr";
        then
            chrNames=$(cut -f1 {input.ucsc2ncbi} | tr '\n' '|')
        else
            chrNames=$(cut -f2 {input.ucsc2ncbi} | tr '\n' '|')
        fi

        # write coverage from idxstats into file
        count=$({params.samtools} idxstats {input.bam} | grep -E "^($chrNames)" | \
                cut -f3 | paste -sd+ - | bc)

        echo -e "{wildcards.sampleID}\t${{count}}" > {output}
        """

rule aberrantExpression_mergeBamStats:
    input: 
        lambda w: expand(cfg.getProcessedDataDir() +
            "/aberrant_expression/{{annotation}}/coverage/{sampleID}.tsv",
            sampleID=sa.getIDsByGroup(w.dataset))
    output:
        cfg.getProcessedDataDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/bam_coverage.tsv"
    params:
        ids = lambda w: sa.getIDsByGroup(w.dataset)
    run:
        with open(output[0], "w") as bam_stats:
            bam_stats.write("sampleID\trecord_count\n")
            for f in input:
                bam_stats.write(open(f, "r").read())

rulegraph_filename = f'{config["htmlOutputPath"]}/AE_rulegraph'

rule aberrantExpression_rulegraph:
    output:
        svg = f"{rulegraph_filename}.svg",
        png = f"{rulegraph_filename}.png"
    shell:
        """
        snakemake aberrantExpression --rulegraph | dot -Tsvg > {output.svg}
        snakemake aberrantExpression --rulegraph | dot -Tpng > {output.png}
        """
