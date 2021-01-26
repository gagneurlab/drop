import pytest
import os
import re

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
            ("splicingCounts", ['Output/processed_results/exported_counts/fraser--hg19--v29/k_j_counts.tsv.gz']),
            ("splicingCounts", ['Output/processed_results/exported_counts/fraser--hg19--v29/k_theta_counts.tsv.gz']),
            ("splicingCounts", ['Output/processed_results/exported_counts/fraser--hg19--v29/n_psi5_counts.tsv.gz']),
            ("splicingCounts", ['Output/processed_results/exported_counts/fraser--hg19--v29/n_theta_counts.tsv.gz'])
        ]
    )
def test_cfgExportCountFiles(demo_dir, dropConfig, prefix, files):
    curFile = files[0]
    true_path = [f"{demo_dir}/{curFile}"]
    if prefix == "geneCounts":
        test_path = dropConfig.exportCounts.getExportCountFiles(prefix)
    else:
        filePattern=os.path.basename(curFile).split(".")[0]
        test_path = dropConfig.exportCounts.getExportCountFiles(prefix, type=[filePattern.split("_")[1]],
                expandPattern=re.sub("_.*_", "_{type}_", filePattern))
    assert true_path == test_path
