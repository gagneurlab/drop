class Test_MAE_Config:

    def test_config(self, demo_dir, dropConfig):
        assert dropConfig.MAE.getWorkdir() == "Scripts/MonoallelicExpression/pipeline"
        dict_ = {
            'qcVcf': f'{demo_dir}/Data/qc_vcf_1000G.vcf.gz',
            'groups': ['mae'],
            'qcGroups': ['mae'],
            'gatkIgnoreHeaderCheck': True,
            'padjCutoff': 0.05,
            'allelicRatioCutoff': 0.8,
            'maxAF': 0.001,
            'addAF': False,
            'maxVarFreqCohort': 0.04,
            'gnomAD': False
        }
        assert dict_.items() <= dropConfig.MAE.dict_.items()

    def test_createMaeIDS(self, dropConfig):
        mae_ids_real = {'mae': ['HG00096--HG00096.1.M_111124_6', 'HG00103--HG00103.4.M_120208_3'].sort()}
        mae_ids_test = {k: l.sort() for k, l in dropConfig.MAE.createMaeIDS().items()}
        assert mae_ids_real == mae_ids_test
