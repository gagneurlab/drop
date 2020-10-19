import pytest


def test_DropConfigKeys(dropConfig):
    for key in dropConfig.CONFIG_KEYS:
        assert dropConfig.get(key) == dropConfig.config_dict[key]


def test_DropConfigPaths(demo_dir, dropConfig):
    assert dropConfig.getRoot() == f"{demo_dir}/Output"
    assert dropConfig.getProcessedDataDir() == f"{demo_dir}/Output/processed_data"
    assert dropConfig.getProcessedResultsDir() == f"{demo_dir}/Output/processed_results"
    gene_anno = {'v29': f'{demo_dir}/Data/gencode_annotation_trunc.gtf'}
    assert gene_anno == dropConfig.getGeneAnnotations()


@pytest.mark.parametrize(
        "modules,groups",
        [
            ("aberrantExpression", {"outrider"}),
            ("aberrantSplicing", {"fraser"}),
            (["aberrantExpression", "aberrantSplicing"], {"outrider", "fraser"})
        ]
    )
def test_cfgExportGroups(dropConfig, modules, groups):
    assert groups == dropConfig.exportCounts.getExportGroups(modules)


@pytest.mark.parametrize(
        "prefix,files",
        [
            ("geneCounts", ['Output/processed_results/exported_counts/outrider--hg19--v29/geneCounts.tsv.gz']),
            ("splitCounts", ['Output/processed_results/exported_counts/fraser--hg19--v29/splitCounts.tsv.gz']),
            ("spliceSiteOverlapCounts", ['Output/processed_results/exported_counts/fraser--hg19--v29/spliceSiteOverlapCounts.tsv.gz'])
        ]
    )
def test_cfgExportCountFiles(demo_dir, dropConfig, prefix, files):
    true_path = [f"{demo_dir}/{file}" for file in files]
    test_path = dropConfig.exportCounts.getExportCountFiles(prefix)
    assert true_path == test_path
