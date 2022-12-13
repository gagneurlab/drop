class Test_AE_Config:

    def test_config(self, dropConfig,demo_dir):
        assert dropConfig.AE.getWorkdir() == f"{demo_dir}/Scripts/AberrantExpression/pipeline"
        dict_ = {
            'groups': ['outrider', 'outrider_external'],
            'fpkmCutoff': 1,
            'implementation': 'autoencoder',
            'padjCutoff': 1,
            'zScoreCutoff': 0,
            'reportAllGenesToTest': False,
            'maxTestedDimensionProportion': 3
        }
        print(dropConfig.AE.dict_.items())
        assert dict_.items() <= dropConfig.AE.dict_.items()

    def test_getCountsFiles(self, demo_dir, dropConfig):
        counts_dir = f"{demo_dir}/Output/processed_data/aberrant_expression/v29/counts"
        counts_files = [
            'HG00103.Rds',
            'HG00096.Rds',
            'HG00106.Rds',
            'HG00116.Rds',
            'HG00111.Rds',
            'HG00149.Rds',
            'HG00176.Rds',
            'HG00150.Rds',
            'HG00126.Rds',
            'HG00132.Rds'
        ]
        counts_files_true = [f"{counts_dir}/{f}" for f in counts_files]
        counts_files_true.sort()
        counts_files_test = dropConfig.AE.getCountFiles(annotation="v29", group="outrider")
        counts_files_test.sort()
        assert counts_files_true == counts_files_test

        # import count
        counts_files_true = counts_files_true[2:]
        counts_files_true.append(f"{demo_dir}/Data/external_count_data/geneCounts.tsv.gz")
        counts_files_true.sort()
        counts_files_test = dropConfig.AE.getCountFiles(annotation="v29", group="outrider_external")
        counts_files_test.sort()
        assert counts_files_true == counts_files_test

    def test_getCountParams(self, dropConfig):
        params = {
            'STRAND': 'no',
            'COUNT_MODE': 'IntersectionStrict',
            'PAIRED_END': True,
            'COUNT_OVERLAPS': True
        }
        print(dropConfig.AE.getCountParams("HG00103"))
        assert params == dropConfig.AE.getCountParams("HG00103")


