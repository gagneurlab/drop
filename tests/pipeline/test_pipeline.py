from tests.common import *


def test_dryrun(demo_dir):
    r = run("snakemake -n --cores 1", dir_path=demo_dir)
    message = "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
    assert message in r.stdout


def test_dependencyGraph(demo_dir):
    r = run(f"snakemake dependencyGraph -F --cores {CORES}", dir_path=demo_dir)
    assert "aberrantExpression_dependency" in r.stderr
    assert "aberrantSplicing_dependency" in r.stderr
    assert "mae_dependency" in r.stderr
    assert "Finished job 0." in r.stderr
    # clean HTML output
    r = run(f"snakemake clean --cores {CORES}", dir_path=demo_dir)
    assert "clean" in r.stderr
    assert "Finished job 0." in r.stderr


def test_export(demo_dir):
    r = run(f"snakemake exportCounts --cores {CORES}", demo_dir)
    assert "Finished job 0." in r.stderr


@pytest.mark.trylast
def test_all(demo_dir):
    r = run(f"snakemake --cores {CORES}", demo_dir)
    assert "Finished job 0." in r.stderr
