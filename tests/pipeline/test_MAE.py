from tests.common import *


class Test_MAE_Pipeline:

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run MAE pipeline")
        pipeline_run = run(["snakemake", "mae", f"-j{CORES}"], demo_dir)
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
