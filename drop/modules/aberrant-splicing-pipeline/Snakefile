WORKDIR = cfg.AS.getWorkdir(str_=False)

rule aberrantSplicing:
    input:
        expand(
            cfg.getHtmlFromScript(WORKDIR / "Counting" / "DatasetsF.R"),
            annotation=cfg.getGeneVersions()
        ),
        expand(
            cfg.getHtmlFromScript(WORKDIR / "FRASER" / "Datasets.R"),
            annotation=cfg.getGeneVersions()
        )

rulegraph_filename = f'{config["htmlOutputPath"]}/AS_rulegraph'
rule aberrantSplicing_rulegraph:
    output:
        svg = f"{rulegraph_filename}.svg",
        png = f"{rulegraph_filename}.png"
    shell:
        """
        snakemake aberrantSplicing --rulegraph | dot -Tsvg > {output.svg}
        snakemake aberrantSplicing --rulegraph | dot -Tpng > {output.png}
        """
