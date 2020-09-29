from tests.conftest import *


class AE_Pipeline:

    def __init__(self, demo_dir):
        self.demo_dir = demo_dir
        print("run aberrant expression pipeline...")
        pipeline_run = run(["snakemake", "aberrantExpression", "-j", CORES], demo_dir)
        assert "Finished job 0." in pipeline_run.stderr

    def counts(self):
        cnt_file = "Output/processed_data/aberrant_expression/v29/outrider/outrider/total_counts.Rds"
        r_cmd = """
            counts <- readRDS("{}")
            print(counts)
            """.format(cnt_file)
        r = runR(r_cmd, self.demo_dir)
        assert "dim: 805 12" in r.stdout
        assert "rownames(805): ENSG00000141956.13_3 ENSG00000141959.16_1 ..." in r.stdout

    def results(self):
        output_dir = "Output/processed_results/aberrant_expression/v29/outrider/outrider"
        r_cmd = """ 
                # ods object
                ods <- readRDS(file.path("{}", "ods.Rds"))
                print(ods)
    
                # results table
                res <- readRDS(file.path("{}", "OUTRIDER_results_all.Rds"))
                print(paste(c("res:", dim(res)), collapse=" "))
                """.format(output_dir, output_dir)
        r = runR(r_cmd, self.demo_dir)
        assert "class: OutriderDataSet" in r.stdout
        assert "dim: 441 12" in r.stdout
        assert "res: 5292 15" in r.stdout