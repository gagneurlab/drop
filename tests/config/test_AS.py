class Test_AS_Config:

    def test_config(self, dropConfig):
        assert dropConfig.AS.getWorkdir() == "Scripts/AberrantSplicing/pipeline"
        dict_ = {
            'groups': ['fraser'],
            'recount': True,
            'longRead': False,
            'keepNonStandardChrs': True,
            'filter': False,
            'minExpressionInOneSample': 20,
            'minDeltaPsi': 0.05,
            'implementation': 'PCA',
            'padjCutoff': 1,
            'zScoreCutoff': 0,
            'deltaPsiCutoff': 0.05,
            'maxTestedDimensionProportion': 6
        }
        assert dict_.items() <= dropConfig.AS.dict_.items()

    def test_getSplitCountFiles(self, demo_dir, dropConfig):
        counts_dir = f"{demo_dir}/Output/processed_data/aberrant_splicing/datasets/cache/raw-fraser/sample_tmp/" \
                     "splitCounts"
        ids = [
            'HG00096.1.M_111124_6_trunc', 'HG00103.4.M_120208_3_trunc', 'HG00111.2.M_111215_4_trunc',
            'HG00106.4.M_120208_5_trunc', 'HG00126.1.M_111124_8_trunc', 'HG00149.1.M_111124_6_trunc',
            'HG00116.2.M_120131_1_trunc', 'HG00132.2.M_111215_4_trunc', 'HG00150.4.M_120208_7_trunc',
            'HG00176.4.M_120208_2_trunc'
        ]
        counts_files_true = [f"{counts_dir}/sample_{id}.done" for id in ids]
        counts_files_test = dropConfig.AS.getSplitCountFiles(dataset="fraser")

        counts_files_true.sort()
        counts_files_test.sort()
        assert counts_files_true == counts_files_test

    def test_getNonSplitCountFiles(self, demo_dir, dropConfig):
        counts_dir = f"{demo_dir}/Output/processed_data/aberrant_splicing/datasets/cache/raw-fraser/sample_tmp/" \
                     "nonSplitCounts"
        ids = [
            'HG00096.1.M_111124_6_trunc', 'HG00103.4.M_120208_3_trunc', 'HG00111.2.M_111215_4_trunc',
            'HG00106.4.M_120208_5_trunc', 'HG00126.1.M_111124_8_trunc', 'HG00149.1.M_111124_6_trunc',
            'HG00116.2.M_120131_1_trunc', 'HG00132.2.M_111215_4_trunc', 'HG00150.4.M_120208_7_trunc',
            'HG00176.4.M_120208_2_trunc'
        ]
        counts_files_true = [f"{counts_dir}/sample_{id}.done" for id in ids]
        counts_files_test = dropConfig.AS.getNonSplitCountFiles(dataset="fraser")

        counts_files_true.sort()
        counts_files_test.sort()
        assert counts_files_true == counts_files_test
