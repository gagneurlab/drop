from tests.conftest import *


class AS_Pipeline:

    def __init__(self, demo_dir):
        self.demo_dir = demo_dir
        print("run aberrant splicing pipeline")
        pipeline_run = run(["snakemake", "aberrantSplicing", "-j", CORES], self.demo_dir)
        assert "Finished job 0." in pipeline_run.stderr

    def counts(self):
        cnt_file = "Output/processed_data/aberrant_splicing/datasets/savedObjects/raw-fraser/fds-object.RDS"
        r_cmd = """ 
            library(FRASER)
            fds <- loadFraserDataSet(file="{}")
            print(fds)
            """.format(cnt_file)
        r = runR(r_cmd, self.demo_dir)
        assert "Number of samples:      10" in r.stdout
        assert "Number of junctions:    81" in r.stdout
        assert "Number of splice sites: 9" in r.stdout

    def results(self):
        results_dir = "Output/processed_data/aberrant_splicing/results"
        r = run(f"wc -l {results_dir}/fraser_results_per_junction.tsv", self.demo_dir)
        assert "88 " in r.stdout
        r = run(f"wc -l {results_dir}/fraser_results.tsv", self.demo_dir)
        assert "1 " in r.stdout
