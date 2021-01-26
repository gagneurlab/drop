from tests.common import *


class Test_AS_Pipeline:

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run aberrant splicing pipeline")
        pipeline_run = run(["snakemake", "aberrantSplicing", f"-j{CORES}"], demo_dir)
        assert "Finished job 0." in pipeline_run.stderr
        return pipeline_run

    @pytest.mark.usefixtures("pipeline_run")
    def test_counts(self, demo_dir):
        cnt_file = "Output/processed_data/aberrant_splicing/datasets/savedObjects/raw-fraser/fds-object.RDS"
        r_cmd = """ 
            library(FRASER)
            fds <- loadFraserDataSet(file="{}")
            print(fds)
            """.format(cnt_file)
        r = runR(r_cmd, demo_dir)
        assert "Number of samples:      10" in r.stdout
        assert "Number of junctions:    81" in r.stdout
        assert "Number of splice sites: 9" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    def test_results(self, demo_dir):
        results_dir = "Output/processed_data/aberrant_splicing/results"
        r = run(f"wc -l {results_dir}/fraser_results_per_junction.tsv", demo_dir)
        assert "87 " in r.stdout
        r = run(f"wc -l {results_dir}/fraser_results.tsv", demo_dir)
        assert "1 " in r.stdout
