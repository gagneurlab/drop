class Test_AS_Config:

    def test_config(self, dropConfig,demo_dir):
        assert dropConfig.AS.getWorkdir() == f"{demo_dir}/Scripts/AberrantSplicing/pipeline"
        dict_ = {
            'run': True,
            'groups': ['fraser', 'fraser_external'],
            'recount': True,
            'longRead': False,
            'keepNonStandardChrs': False,
            'filter': False,
            'minExpressionInOneSample': 20,
            'quantileMinExpressionIn': 10,
            'quantileForFiltering': 0.95,
            'minDeltaPsi': 0.05,
            'implementation': 'PCA',
            'padjCutoff': 1,
            'deltaPsiCutoff': 0.05,
            'maxTestedDimensionProportion': 6,
            'FRASER_version': 'FRASER'
        }
        assert dict_.items() <= dropConfig.AS.dict_.items()

    def test_getSplitCountFiles(self, demo_dir, dropConfig):
        counts_dir = f"{demo_dir}/Output/processed_data/aberrant_splicing/datasets/cache/raw-local-fraser/sample_tmp/" \
                     "splitCounts"
        ids = [
            'HG00096', 'HG00103', 'HG00111',
            'HG00106', 'HG00126', 'HG00149',
            'HG00116', 'HG00132', 'HG00150',
            'HG00176'
        ]
        counts_files_true = [f"{counts_dir}/sample_{id}.done" for id in ids]
        counts_files_test = dropConfig.AS.getSplitCountFiles(dataset="fraser")

        counts_files_true.sort()
        counts_files_test.sort()
        assert counts_files_true == counts_files_test

    def test_getNonSplitCountFiles(self, demo_dir, dropConfig):
        counts_dir = f"{demo_dir}/Output/processed_data/aberrant_splicing/datasets/cache/raw-local-fraser/sample_tmp/" \
                     "nonSplitCounts"
        ids = [
            'HG00096', 'HG00103', 'HG00111',
            'HG00106', 'HG00126', 'HG00149',
            'HG00116', 'HG00132', 'HG00150',
            'HG00176'
        ]
        counts_files_true = [f"{counts_dir}/sample_{id}.done" for id in ids]
        counts_files_test = dropConfig.AS.getNonSplitCountFiles(dataset="fraser")

        counts_files_true.sort()
        counts_files_test.sort()
        assert counts_files_true == counts_files_test
