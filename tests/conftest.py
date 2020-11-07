import os
import wbuild
import drop
import sys
from unittest.mock import patch
from .common import *


@pytest.yield_fixture(scope="session", autouse=True)
def demo_dir(tmpdir_factory):
    """
    Directory containing files downloaded from Gagneurlab ready to run the demo pipeline
    :param testdirectory: inherits from pytest-testdirectory
    :return: demo directory
    """
    run_dir = tmpdir_factory.mktemp("demo_dir")
    LOGGER.info(f"\n create demo dir: {run_dir}")
    r = run(["drop", "demo"], run_dir, stdout=subprocess.DEVNULL)
    assert "demo project created" in r.stderr
    yield run_dir
    LOGGER.info("\n remove demo directory")
    run_dir.remove()


@pytest.fixture(scope="session", autouse=True)
def dropConfig(demo_dir):
    orig_path = os.getcwd()
    os.chdir(demo_dir)
    with patch.object(sys, 'argv', ["snakemake", "-n"]):
        cfg = drop.config.DropConfig(wbuild.utils.Config())
    os.chdir(orig_path)
    return cfg
