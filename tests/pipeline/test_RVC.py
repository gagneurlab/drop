from tests.common import *


class Test_RVC_Pipeline:

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
                vcf  <- read.table("{}", stringsAsFactors = FALSE) 
                print(nrow(vcf))
                """.format(vcf_file)
        r = runR(r_cmd, demo_dir)
        assert "[1] 4965" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    # count the number of variant calls in batch0 before splitting
    def test_variants_single_line_multi(self, demo_dir):
        vcf_file = "Output/processed_data/rnaVariantCalling/out/all_samples_haplocaller/batch_0_all_samples.genotyped.filtered_clean.vcf.gz"
        r_cmd = """ 
                library(data.table)
                vcf  <- read.table("{}", stringsAsFactors = FALSE) 
                print(nrow(vcf))
                """.format(vcf_file)
        r = runR(r_cmd, demo_dir)
        assert "[1] 5005" in r.stdout

    @pytest.mark.usefixtures("pipeline_run")
    def test_single_sample_variants(self, demo_dir):
        result_dir = "Output/processed_data/rnaVariantCalling/out/sample_haplocaller/HG00096.1.M_111124_6"
        r_cmd = """ 
                library(data.table)
                basic_filter <- fread(file.path("{}", "HG00096.1.M_111124_6.genotyped.filtered.basic10.vcf.gz"))
                masked_filter <- fread(file.path("{}", "HG00096.1.M_111124_6.genotyped.filtered.basic10.masked.vcf.gz"))

                nrow_basic <- nrow(basic_filter)
                nrow_masked <- nrow(masked_filter)

                print(nrow_basic)
                print(nrow_basic == nrow_masked)
                print(sum(basic_filter[,7] == "PASS"))
                print(sum(masked_filter[,7] == "PASS"))
                """.format(result_dir,result_dir)
        r = runR(r_cmd, demo_dir)
        assert "[1] 3006" in r.stdout
        assert "[1] TRUE" in r.stdout
        assert "[1] 515" in r.stdout
        assert "[1] 482" in r.stdout
