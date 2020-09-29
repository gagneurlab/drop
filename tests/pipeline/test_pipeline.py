from .aberrantExpression import *
from .aberrantSplicing import *
from .mae import *


def test_dryrun(demo_dir):
    r = run(["snakemake", "-n"], dir_path=demo_dir)
    message = "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
    assert message in r.stdout


def test_dependencyGraph(demo_dir):
    r = run(["snakemake", "dependencyGraph", "-j"], dir_path=demo_dir)
    assert "aberrantExpression_dependency" in r.stderr
    assert "aberrantSplicing_dependency" in r.stderr
    assert "mae_dependency" in r.stderr
    assert "Finished job 0." in r.stderr
    # clean HTML output
    r = run(["snakemake", "clean", "-j"], dir_path=demo_dir)
    assert "clean" in r.stderr
    assert "Finished job 0." in r.stderr


def test_export(demo_dir):
    r = run(["snakemake", "exportCounts", "-j", CORES], demo_dir)
    assert "Finished job 0." in r.stderr


def test_submodules(demo_dir):
    test_AE = AE_Pipeline(demo_dir)
    test_AE.counts()
    test_AE.results()

    test_AS = AS_Pipeline(demo_dir)
    test_AS.counts()
    test_AS.results()

    test_MAE = MAE_Pipeline(demo_dir)
    test_MAE.counts()
    test_MAE.results()


def test_all(demo_dir):
    r = run(["snakemake", "-j", CORES], demo_dir)
    assert "Finished job 0." in r.stderr
