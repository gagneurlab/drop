from tests.common import *


class Test_RVC_Pipeline:

# change the third instance of "run: true" to "run: false" to turn off the RVC module
# run the RVC module using this config_RVC_norun (which should do nothing)
    def test_pipeline_no_run(self,demo_dir):
        LOGGER.info("run RNA variant calling pipeline with \'run: false\'")
        # change the fourth instance of "run: true" to "run: false" to turn off the RVC module
        # run the RVC module using this config_RVC_norun (which should do nothing)
        run("awk -v n=4 \'/run: true/ { if (++count == n) sub(/run: true/, \"run: false\"); } 1\' \
        	config.yaml > config_RVC_norun.yaml  ",demo_dir)
        try:
            pipeline_run = run(["snakemake", "rnaVariantCalling", f"-np", "--configfile", "config_RVC_norun.yaml"], demo_dir)
            pipeline_run = run(["snakemake", "rnaVariantCalling", f"-c{CORES}", "--configfile", "config_RVC_norun.yaml"], demo_dir)
        except subprocess.CalledProcessError:
            print("Failed Successfully")

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run RVC pipeline")
        pipeline_run = run(f"snakemake rnaVariantCalling -c{CORES}", demo_dir)
        assert "Finished job 0." in pipeline_run.stderr
        return pipeline_run

    @pytest.mark.usefixtures("pipeline_run")
    # count the number of variant calls in batch0 before splitting
    def test_variants_batch(self, demo_dir):
        vcf_file0 = f"{demo_dir}/Output/processed_results/rnaVariantCalling/out/batch_vcfs/batch_0/batch_0_v29.annotated.vcf.gz"
        vcf_file1 = f"{demo_dir}/Output/processed_results/rnaVariantCalling/out/batch_vcfs/batch_1/batch_1_v29.annotated.vcf.gz"
        r_cmd = """ 
                library(data.table)
                vcf0  <- fread("{}")
                num_variants0 <- nrow(vcf0)
                print(num_variants0)
                vcf1  <- fread("{}")
                num_variants1 <- nrow(vcf1)
                print(num_variants1)
                """.format(vcf_file0,vcf_file1)
        r = runR(r_cmd, demo_dir)
        assert "[1] 585" in r.stdout
        assert "[1] 812" in r.stdout
