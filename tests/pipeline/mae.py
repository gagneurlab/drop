from tests.conftest import *


class MAE_Pipeline:

    def __init__(self, demo_dir):
        self.demo_dir = demo_dir
        print("run MAE pipeline")
        pipeline_run = run(["snakemake", "mae", "-j", CORES], demo_dir)
        assert "Finished job 0." in pipeline_run.stderr

    def counts(self):
        cnt_file = "Output/processed_data/mae/allelic_counts/HG00103--HG00103.4.M_120208_3.csv.gz"
        r_cmd = """ 
                library(data.table)
                cnts <- fread("{}")
                print(nrow(cnts))
                """.format(cnt_file)
        r = runR(r_cmd, self.demo_dir)
        assert "[1] 235" in r.stdout

    def results(self):
        results_file = "Output/processed_results/mae/mae/MAE_results_all_v29.tsv.gz"
        r_cmd = """ 
                library(data.table)
                res <- fread("{}")
                print(nrow(res))
                """.format(results_file)
        r = runR(r_cmd, self.demo_dir)
        assert "335 " in r.stdout
