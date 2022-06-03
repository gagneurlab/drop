class Test_RVC_Config:
    def test_config(self, demo_dir,dropConfig):
        assert dropConfig.RVC.getWorkdir() == demo_dir + "/Scripts/rnaVariantCalling/pipeline"
        dict_ = {
            'groups': ['batch_0','batch_1'],
            'highQualityVCFs': [f'{demo_dir}/Data/high_confidence_snps.vcf.gz', f'{demo_dir}/Data/high_confidence_indels.vcf.gz'],
            'dbSNP': f'{demo_dir}/Data/dbSNP_chr21.vcf.gz',
            'repeat_mask': f'{demo_dir}/Data/repeat_mask_chr21.bed',
            'hcArgs': '',
            'createSingleVCF': True,
            'minAlt': 3
        }
        assert dict_.items() <= dropConfig.RVC.dict_.items()

