class Test_MAE_Config:

    def test_config(self,dropConfig,demo_dir):
        assert dropConfig.MAE.getWorkdir() == f"{demo_dir}/Scripts/MonoallelicExpression/pipeline"
        dict_ = {
            'qcVcf': f'{demo_dir}/Data/qc_vcf_1000G.vcf.gz',
            'groups': ['mae'],
            'qcGroups': ['mae'],
            'gatkIgnoreHeaderCheck': True,
            'padjCutoff': 0.5,
            'allelicRatioCutoff': 0.7,
            'maxAF': 0.001,
            'addAF': False,
            'maxVarFreqCohort': 1,
        }
        assert dict_.items() <= dropConfig.MAE.dict_.items()

    def test_createMaeIDS(self, dropConfig):
        mae_ids_real = {'mae': ['HG00096--HG00096', 'HG00103--HG00103'].sort()}
        mae_ids_test = {k: l.sort() for k, l in dropConfig.MAE.createMaeIDS().items()}
        assert mae_ids_real == mae_ids_test
