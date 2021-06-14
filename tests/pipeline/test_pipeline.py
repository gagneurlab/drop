from tests.common import *


def test_dryrun(demo_dir):
    r = run(["snakemake", "-n"], dir_path=demo_dir)
    message = "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
    assert message in r.stdout

def test_pipeline_no_run(self,demo_dir):
    LOGGER.info("run entire pipeline with \'run: false\'")
    run("sed 's/run: true/run: false:/g' config.yaml > config_norun.yaml  ",demo_dir)
    pipeline_run = run(["snakemake",  f"-j{CORES}", "--configfile", "config_norun.yaml"], demo_dir)
    tmp = run(["snakemake", "--unlock"], demo_dir)
    assert "Finished job 0." in pipeline_run.stderr
    return pipeline_run


def test_dependencyGraph(demo_dir):
    r = run(["snakemake", "dependencyGraph", "-Fj"], dir_path=demo_dir)
    assert "aberrantExpression_dependency" in r.stderr
    assert "aberrantSplicing_dependency" in r.stderr
    assert "mae_dependency" in r.stderr
    assert "Finished job 0." in r.stderr
    # clean HTML output
    r = run(["snakemake", "clean", f"-j{CORES}"], dir_path=demo_dir)
    assert "clean" in r.stderr
    assert "Finished job 0." in r.stderr


def test_export(demo_dir):
    r = run(["snakemake", "exportCounts", f"-j{CORES}"], demo_dir)
    assert "Finished job 0." in r.stderr


@pytest.mark.trylast
def test_all(demo_dir):
    r = run(["snakemake", f"-j{CORES}"], demo_dir)
    assert "Finished job 0." in r.stderr
