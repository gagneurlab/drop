import drop
from drop.config.DropConfig import DropConfig
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


def installRPackages(config: DropConfig = None):
    logger.info("check for missing R packages")
    script = Path(drop.__file__).parent / "installRPackages.R"
    requirements = Path(drop.__file__).parent / 'requirementsR.txt'

    # install main packages
    response = subprocess.run(["Rscript", script, requirements], stderr=subprocess.STDOUT)
    response.check_returncode()

    # install pipeline depending packages
    if config is not None:
        pkg_assembly_name = None
        pkg_mafdb_name = None

        if config.genomeAssembly in ['hg19']:
            pkg_assembly_name = "BSgenome.Hsapiens.UCSC.hg19"
            if config.mae.addAF:
                pkg_mafdb_name = "MafDb.gnomAD.r2.1.hs37d5"
        else if config.genomeAssembly in ['hg38']
            pkg_assembly_name = "BSgenome.Hsapiens.UCSC.hg38"
            if config.mae.addAF:
                pkg_mafdb_name = "MafDb.gnomAD.r2.1.GRCh38"

        if pkg_name is not None:
            response = subprocess.run(["Rscript", script, pkg_name], stderr=subprocess.STDOUT)
            response.check_returncode()

        if pkg_mafdb_name is not None:
            response = subprocess.run(["Rscript", script, pkg_mafdb_name], stderr=subprocess.STDOUT)
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
