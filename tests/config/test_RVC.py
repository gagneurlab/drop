class Test_RVC_Config:
    def test_config(self, demo_dir,dropConfig):
        assert dropConfig.RVC.getWorkdir() == demo_dir + "/Scripts/rnaVariantCalling/pipeline"
        dict_ = {
            'groups': ['batch_0','batch_1'],
            'knownVCFs': [f'{demo_dir}/Data/high_confidence_snps_test.vcf.gz', f'{demo_dir}/Data/high_confidence_indels_test.vcf.gz', f'{demo_dir}/Data/dbSNP_chr21_test.vcf.gz'],
            'repeat_mask': f'{demo_dir}/Data/repeat_mask_chr21_test.bed',
            'hcArgs': '',
            'minAlt': 10
        }
        print(dropConfig.RVC.dict_.items())
        assert dict_.items() <= dropConfig.RVC.dict_.items()

