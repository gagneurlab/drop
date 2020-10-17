from click.testing import CliRunner
from drop import cli
from .common import *


def test_Help_isShown():
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help' in help_result.output
    assert 'Show this message and exit.' in help_result.output


def test_drop_init(tmpdir):
    init_dir = tmpdir.mkdir("init")
    r = run(["drop", "init"], dir_path=init_dir)
    assert 'init...done\n' in r.stderr


# TOOD: test update
