import drop
import subprocess
from pathlib import Path
from snakemake.logging import logger


def setupPaths(projectRoot):
    # repository paths
    repoRoot = Path(drop.__file__).parent
    repoPaths = {
        "template": repoRoot / "template",
        "Scripts": repoRoot / "template" / "Scripts",
        "modules": repoRoot / "modules"
    }

    # project paths
    projectPaths = {
        "projectDir": projectRoot,
        "Scripts": projectRoot / "Scripts",
        "dropDir": projectRoot / ".drop",
        "tmpDir": projectRoot / ".drop" / "tmp"
    }

    return repoPaths, projectPaths


def installRPackages():
    logger.info("check for missing R packages")
    script = Path(drop.__file__).parent / "installRPackages.R"
    requirements = Path(drop.__file__).parent / 'requirementsR.txt'

    response = subprocess.run(["Rscript", script, requirements], stderr=subprocess.STDOUT)
    response.check_returncode()


def checkDropVersion(projectRoot, force=False):
    if projectRoot != Path.cwd().resolve():
        raise AssertionError(f"Specified project root '{projectRoot}' not equal to current working directory "
                             f"'{Path.cwd().resolve()}'")

    version_file = projectRoot / ".drop"/ ".version"
    version_file.touch()
    with open(version_file, "r") as f:
        project_version = f.readline()

    if drop.__version__ != project_version:
        logger.info(f"Update drop version for {projectRoot} to version {drop.__version__}")
        with open(version_file, "w") as f:
            f.write(drop.__version__)
        drop.cli.setFiles(projectRoot)
    elif force:
        drop.cli.setFiles(projectRoot)
