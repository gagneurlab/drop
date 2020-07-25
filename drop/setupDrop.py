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
        "modules": projectRoot / ".drop" / "modules",
        "tmpDir": projectRoot / ".drop" / "tmp"
    }

    return repoPaths, projectPaths

def installRPackages():
    logger.info("check for missing R packages")
    script = Path(drop.__file__).parent / "installRPackages.R"
    requirements = Path(drop.__file__).parent / 'requirementsR.txt'
    
    response = subprocess.run(["Rscript", script, requirements], stderr=subprocess.STDOUT)
    response.check_returncode()
