from click.testing import CliRunner
from drop import cli
from .conftest import *


cores = 2


def test_Help_isShown():
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help' in help_result.output
    print(help_result.output)
    assert 'Show this message and exit.' in help_result.output


def test_drop_init(tmpdir):
    init_dir = tmpdir.mkdir("init")
    r = run(["drop", "init"], dir_path=init_dir)
    assert 'init...done\n' == r.stderr


def test_dryrun(demo_dir):
    r = run(["snakemake", "-n"], dir_path=demo_dir)
    message = "This was a dry-run (flag -n). The order of jobs does not reflect the order of execution."
    assert message in r.stdout


def test_aberrantExpression(demo_dir):
    r = run(["snakemake", "aberrantExpression", f"-j{cores}"], demo_dir)
    assert "Finished job 0." in r.stderr


def test_aberrantSplicing(demo_dir):
    r = run(["snakemake", "aberrantSplicing", f"-j{cores}"], demo_dir)
    assert "Finished job 0." in r.stderr

def test_MAE(demo_dir):
    r = run(["snakemake", "mae", f"-j{cores}"], demo_dir)
    assert "Finished job 0." in r.stderr

