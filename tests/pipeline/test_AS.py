from tests.common import *


class Test_AS_Pipeline:

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run aberrant splicing pipeline")
        pipeline_run = run(["snakemake", "aberrantSplicing", f"-j{CORES}"], demo_dir)
        tmp = run(["snakemake", "--unlock"], demo_dir)
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
        assert "Number of junctions:    81" in r.stdout
        assert "Number of splice sites: 9" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    def test_results(self, demo_dir):
        results_dir = "Output/processed_results/aberrant_splicing/results"                                              
        annotation = "v29"                                                                                              
        dataset = "fraser"                                                                                              
        r = run(f"wc -l {results_dir}/{annotation}/fraser/{dataset}/results_per_junction.tsv", demo_dir)                
        assert "87" == r.stdout.split()[0]                                                                                        
        r = run(f"wc -l {results_dir}/{annotation}/fraser/{dataset}/results.tsv", demo_dir)                             
        assert "11" == r.stdout.split()[0]                                     
