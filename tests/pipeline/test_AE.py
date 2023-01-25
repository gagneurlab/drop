from tests.common import *


class Test_AE_Pipeline:

    def test_pipeline_no_run(self,demo_dir):
        flag = False
        LOGGER.info("run aberrantExpression pipeline with \'run: false\'")
        # change the first instance of "run: true" to "run: false" to turn off the AE module
        # run the AE module using this config_AE_norun (which should do nothing)
        run("awk -v n=1 \'/run: true/ { if (++count == n) sub(/run: true/, \"run: false\"); } 1\' \
          config.yaml > config_AE_norun.yaml  ",demo_dir)
        try:
            pipeline_run = run(["snakemake", "aberrantExpression", f"-c{CORES}", "--configfile", "config_AE_norun.yaml"], demo_dir)
        except subprocess.CalledProcessError:
            print("Failed Successfully")

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run aberrant expression pipeline...")
        pipeline_run = run(f"snakemake aberrantExpression -c{CORES}", demo_dir)
        assert "Finished job 0." in pipeline_run.stderr
        return pipeline_run

    @pytest.mark.usefixtures("pipeline_run")
    def test_counts(self, demo_dir):
        cnt_file = "Output/processed_data/aberrant_expression/v29/outrider/outrider/total_counts.Rds"
        r_cmd = """
            counts <- readRDS("{}")
            print(counts)
            """.format(cnt_file)
        r = runR(r_cmd, demo_dir)
        assert "dim: 805 10" in r.stdout
        assert "rownames(805): ENSG00000141956.13_3 ENSG00000141959.16_1 ..." in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    def test_results(self, demo_dir):
        output_dir = "Output/processed_results/aberrant_expression/v29/outrider/outrider"
        r_cmd = """
                # ods object
                ods <- readRDS(file.path("{}", "ods.Rds"))
                print(ods)

                # results table
                res <- readRDS(file.path("{}", "OUTRIDER_results_all.Rds"))
                print(paste(c("res:", dim(res)), collapse=" "))
                """.format(output_dir, output_dir)
        r = runR(r_cmd, demo_dir)
        assert "class: OutriderDataSet" in r.stdout
        assert "dim: 161 10" in r.stdout
        assert "res: 1820 17" in r.stdout

    def test_import_results(self, demo_dir):
        output_dir = "Output/processed_results/aberrant_expression/v29/outrider/outrider_external"
        r_cmd = """
                # ods object
                ods <- readRDS(file.path("{}", "ods.Rds"))
                print(ods)

                # results table
                res <- readRDS(file.path("{}", "OUTRIDER_results_all.Rds"))
                print(paste(c("res:", dim(res)), collapse=" "))
                """.format(output_dir, output_dir)
        r = runR(r_cmd, demo_dir)
        assert "class: OutriderDataSet" in r.stdout
        assert "dim: 389 10" in r.stdout
        assert "res: 4050 17" in r.stdout

    @pytest.fixture()
    def no_import(self, demo_dir):
        LOGGER.info("dryrun without import counts...")

        # adapt config
        run("sed '/outrider_external/d' config.yaml > config_noimp.yaml", demo_dir)

        yield demo_dir

        # reset changed files back to original
        run("rm config_noimp.yaml", demo_dir)

    def test_no_import(self, no_import):
        merged_counts = f"{no_import}/Output/processed_data/aberrant_expression/v29/outrider/outrider/total_counts.Rds"
        # check dryrun
        r = run(f"snakemake {merged_counts} --configfile config_noimp.yaml -nF -c1", no_import)

        # check if new rules are included in snakemake output
        assert "AberrantExpression_pipeline_Counting_countReads_R" in r.stdout
        assert "AberrantExpression_pipeline_Counting_mergeCounts_R" in r.stdout

        # check if pipeline runs through without errors
        run(f"snakemake {merged_counts} --configfile config_noimp.yaml -c{CORES}", no_import)
