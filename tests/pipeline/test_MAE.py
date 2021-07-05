from tests.common import *


class Test_MAE_Pipeline:

    def test_pipeline_no_run(self,demo_dir):
        LOGGER.info("run monoallelicExpression pipeline with \'run: false\'")
        # change the third instance of "run: true" to "run: false" to turn off the MAE module
        # run the MAE module using this config_MAE_norun (which should do nothing)
        run("awk -v n=3 \'/run: true/ { if (++count == n) sub(/run: true/, \"run: false\"); } 1\' \
          config.yaml > config_MAE_norun.yaml  ",demo_dir)
        pipeline_run = run(["snakemake", "aberrantExpression", f"-j{CORES}", "--configfile", "config_MAE_norun.yaml"], demo_dir)
        tmp = run(["snakemake", "--unlock"], demo_dir)
        assert "Nothing to be done." in pipeline_run.stderr
        return pipeline_run

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run MAE pipeline")
        pipeline_run = run(f"snakemake mae --cores {CORES}", demo_dir)
        assert "Finished job 0." in pipeline_run.stderr
        return pipeline_run

    @pytest.mark.usefixtures("pipeline_run")
    def test_counts(self, demo_dir):
        cnt_file = "Output/processed_data/mae/allelic_counts/HG00103--HG00103.4.M_120208_3.csv.gz"
        r_cmd = """
                library(data.table)
                cnts <- fread("{}")
                print(nrow(cnts))
                """.format(cnt_file)
        r = runR(r_cmd, demo_dir)
        assert "[1] 235" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    def test_results(self, demo_dir):
        results_file = "Output/processed_results/mae/mae/MAE_results_all_v29.tsv.gz"
        r_cmd = """
                library(data.table)
                res <- fread("{}")
                print(nrow(res))
                """.format(results_file)
        r = runR(r_cmd, demo_dir)
        assert "[1] 335" in r.stdout
