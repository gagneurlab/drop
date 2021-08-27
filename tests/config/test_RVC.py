class Test_RVC_Config:
    def test_config(self, demo_dir,dropConfig):
        assert dropConfig.RVC.getWorkdir() == demo_dir + "/Scripts/rnaVariantCalling/pipeline"
        dict_ = {
            'groups': ['batch_0','batch_1'],
            'knownVCFs': ['Data/high_confidence_snps_test.vcf.gz', 'Data/high_confidence_indels_test.vcf.gz', 'Data/dbSNP_chr21_test.vcf.gz'],
            'repeat_mask': 'Data/repeat_mask_chr21_test.bed',
            'hcArgs': '',
            'minAlt': 10
        }
        print(dropConfig.RVC.dict_.items())
        assert dict_.items() <= dropConfig.RVC.dict_.items()

