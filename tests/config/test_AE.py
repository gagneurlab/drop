class Test_AE_Config:

    def test_config(self, dropConfig):
        assert dropConfig.AE.getWorkdir() == "Scripts/AberrantExpression/pipeline"
        dict_ = {
            'groups': ['outrider', 'import_exp'],
            'fpkmCutoff': 1,
            'implementation': 'autoencoder',
            'padjCutoff': 1,
            'zScoreCutoff': 0,
            'maxTestedDimensionProportion': 3
        }
        print(dropConfig.AE.dict_.items())
        assert dict_.items() <= dropConfig.AE.dict_.items()

    def test_getCountsFiles(self, demo_dir, dropConfig):
        counts_dir = f"{demo_dir}/Output/processed_data/aberrant_expression/v29/counts"
        counts_files = [
            'HG00103.4.M_120208_3.Rds',
            'HG00096.1.M_111124_6.Rds',
            'HG00106.4.M_120208_5.Rds',
            'HG00116.2.M_120131_1.Rds',
            'HG00111.2.M_111215_4.Rds',
            'HG00149.1.M_111124_6.Rds',
            'HG00176.4.M_120208_2.Rds',
            'HG00150.4.M_120208_7.Rds',
            'HG00126.1.M_111124_8.Rds',
            'HG00132.2.M_111215_4.Rds'
        ]
        counts_files_true = [f"{counts_dir}/{f}" for f in counts_files]
        counts_files_true.sort()
        counts_files_test = dropConfig.AE.getCountFiles(annotation="v29", group="outrider")
        counts_files_test.sort()
        assert counts_files_true == counts_files_test

        # import count
        counts_files_true = counts_files_true[2:]
        counts_files_true.append(f"{demo_dir}/Data/external_geneCounts.tsv.gz")
        counts_files_true.sort()
        counts_files_test = dropConfig.AE.getCountFiles(annotation="v29", group="import_exp")
        counts_files_test.sort()
        assert counts_files_true == counts_files_test

    def test_getCountParams(self, dropConfig):
        params = {
            'STRAND': 'no',
            'COUNT_MODE': 'IntersectionStrict',
            'PAIRED_END': True,
            'COUNT_OVERLAPS': True
        }
        assert params == dropConfig.AE.getCountParams("HG00103.4.M_120208_3")


