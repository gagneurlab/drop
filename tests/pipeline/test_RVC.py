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
            pipeline_run = run(["snakemake", "aberrantExpression", f"-j{CORES}", "--configfile", "config_RVC_norun.yaml"], demo_dir)
            raise AssertionError
        except subprocess.CalledProcessError:
            print("Failed Successfully")

    @pytest.fixture(scope="class")
    def pipeline_run(self, demo_dir):
        LOGGER.info("run RVC pipeline")
        pipeline_run = run(["snakemake", "rnaVariantCalling", f"-j{CORES}"], demo_dir)
        assert "Finished job 0." in pipeline_run.stderr
        return pipeline_run

    @pytest.mark.usefixtures("pipeline_run")
    # count the number of variant calls in batch0 before splitting
    def test_variants_joint(self, demo_dir):
        vcf_file = "Output/processed_data/rnaVariantCalling/out/all_samples_haplocaller/batch_0_all_samples.genotyped.vcf.gz"
        r_cmd = """ 
                library(data.table)
                vcf  <- fread("{}")
                num_variants <- nrow(vcf)
                print(num_variants)
                """.format(vcf_file)
        r = runR(r_cmd, demo_dir)
        assert "[1] 4966" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    # count the number of variant calls in batch0 before splitting
    def test_variants_single_line_multi(self, demo_dir):
        vcf_file = "Output/processed_data/rnaVariantCalling/out/all_samples_haplocaller/batch_0_all_samples.genotyped.filtered_clean.vcf.gz"
        r_cmd = """ 
                library(data.table)
                vcf  <- fread("{}")
                num_variants <- nrow(vcf)
                print(num_variants)
                """.format(vcf_file)
        r = runR(r_cmd, demo_dir)
        assert "[1] 5006" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    def test_single_sample_variants(self, demo_dir):
        result_dir = "Output/processed_data/rnaVariantCalling/out/sample_haplocaller/HG00096"
        r_cmd = """ 
                library(data.table)
                masked_filter <- fread(file.path("{}", "HG00096.genotyped.filtered.basic10.masked.vcf.gz"))

                nrow_masked <- nrow(masked_filter)

                print(sum(masked_filter[,7] == "PASS"))
                """.format(result_dir,result_dir)
        r = runR(r_cmd, demo_dir)
        assert "[1] 482" in r.stdout
