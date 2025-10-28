from tests.common import *


class Test_AS_Pipeline:

    def test_pipeline_no_run(self,demo_dir):
        LOGGER.info("run aberrantSplicing pipeline with \'run: false\'")
        # change the second instance of "run: true" to "run: false" to turn off the AS module
        # run the AS module using this config_AS_norun (which should do nothing)
        run("awk -v n=2 \'/run: true/ { if (++count == n) sub(/run: true/, \"run: false\"); } 1\' \
          config.yaml > config_AS_norun.yaml  ",demo_dir)
        try:
            pipeline_run = run(["snakemake", "aberrantSplicing", f"-c{CORES}", "--configfile", "config_AS_norun.yaml"], demo_dir)
        except subprocess.CalledProcessError:
            print("Failed Successfully")

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run aberrant splicing pipeline")
        pipeline_run = run(f"snakemake aberrantSplicing -c{CORES}", demo_dir)

        assert "Finished job 0." in pipeline_run.stderr
        return pipeline_run


#once Output is prebuilt with ncbi_fds obj present
#    @pytest.mark.usefixtures("pipeline_run")
#    def pipeline_run(self, demo_dir):
#        LOGGER.info("run aberrant splicing ncbi results")
#        annotation = "v29"
#        dataset = "fraser_ncbi"
#        pipeline_run = run(["snakemake", f"{demo_dir}/Output/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv",
#                            f"-j{CORES}"], demo_dir)
#        message = "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
#        assert message in r.stdout


    @pytest.mark.usefixtures("pipeline_run")
    def test_counts(self, demo_dir):
        annotation = "v29"
        dataset = "fraser"
        cnt_file = f"Output/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/fds-object.RDS"
        r_cmd = """
            library(FRASER)
            fds <- loadFraserDataSet(file="{}")
            print(fds)
            """.format(cnt_file)
        r = runR(r_cmd, demo_dir)
        assert "Number of samples:      10" in r.stdout
        assert "Number of junctions:    982" in r.stdout
        assert "Number of splice sites: 356" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    def test_results(self, demo_dir):
        results_dir = "Output/processed_results/aberrant_splicing/results"
        annotation = "v29"
        dataset = "fraser"
        r = run(f"wc -l {results_dir}/{annotation}/fraser/{dataset}/results_per_junction.tsv", demo_dir)
        assert "848" == r.stdout.split()[0] 
        r = run(f"wc -l {results_dir}/{annotation}/fraser/{dataset}/results.tsv", demo_dir)
        assert "259" == r.stdout.split()[0]
        r = run(f"wc -l {results_dir}/{annotation}/fraser/{dataset}/results_gene_all.tsv", demo_dir)
        assert "1701" == r.stdout.split()[0]

