import pytest
import subprocess
import os
import logging

LOGGER = logging.getLogger(__name__)

@pytest.yield_fixture(scope="session")
def demo_dir(tmpdir_factory):
    """
    Directory containing files downloaded from Gagneurlab ready to run the demo pipeline
    :param testdirectory: inherits from pytest-testdirectory
    :return: demo directory
    """
    run_dir = tmpdir_factory.mktemp("demo_dir")
    print(f"\n create demo dir: {run_dir}")
    r = run(["drop", "demo"], run_dir)
    assert 'demo project created' in r.stderr
    yield run_dir
    print("\n remove demo directory")
    run_dir.remove()


def run(cmd, dir_path, report_stdout=True, **kwargs):
    shell = isinstance(cmd, str)
    response = subprocess.run(
        cmd, cwd=dir_path, shell=shell,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        **kwargs
    )
    try:
        response.check_returncode()
    except subprocess.CalledProcessError:
        if report_stdout:
            LOGGER.error("Standard out:")
            LOGGER.error(response.stdout)
        LOGGER.error(response.stderr)
        raise
    return response


def runR(r_cmd, dir_path, report_stdout=False):
    return run(f"Rscript -e '{r_cmd}'", dir_path, report_stdout=report_stdout)

"""
def run(cmd, dir_path, **kwargs):
    shell = isinstance(cmd, str)
    os.environ['PYTHONUNBUFFERED'] = "1"
    proc = subprocess.Popen(
        cmd,
        cwd=dir_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=False,
        **kwargs
    )
    stdout = []
    stderr = []
    mix = []
    while proc.poll() is None:
        line = proc.stdout.readline()
        if line != "":
            stdout.append(line)
            mix.append(line)
            print(line, end='')

        line = proc.stderr.readline()
        if line != "":
            stderr.append(line)
            mix.append(line)
            print(line, end='')

    return proc.returncode, stdout, stderr, mix
"""