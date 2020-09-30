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
