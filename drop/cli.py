import wbuild
import drop
import yaml
from pathlib import Path
from shutil import copyfile
from distutils.dir_util import mkpath, copy_tree, remove_tree
import subprocess
import click
import click_log
import logging

wbuildPath = Path(wbuild.__file__).parent / ".wBuild"
logger = logging.getLogger(__name__)
click_log.basic_config(logger)

@click.group()
@click_log.simple_verbosity_option(logger)
@click.version_option('1.0.2',prog_name='drop')
def main():
    pass


def copyModuleCode(repoPaths, projectPaths):
    repo_map = {
        "aberrant-expression-pipeline": "AberrantExpression",
        "aberrant-splicing-pipeline": "AberrantSplicing",
        "mae-pipeline": "MonoallelicExpression",
        "rvc-pipeline": "rnaVariantCalling"
    }

    for repo, analysis_dir in repo_map.items():
        module_repo = repoPaths["modules"] / repo
        module_project = projectPaths["Scripts"] / analysis_dir / "pipeline"
        if module_project.is_dir():
            remove_tree(module_project)
            print(f"overwriting pipeline scripts for {analysis_dir}")
        copy_tree(str(module_repo), str(module_project))


def removeFile(filePath, warn=True):
    filePath = Path(filePath)
    if filePath.is_file():
        if warn:
            input(f"The file {str(filePath)} will be overwritten. Press Enter to continue...")
        filePath.unlink()


def setFiles(projectDir=None, warn=True):
    projectDir = Path.cwd().resolve() if projectDir is None else projectDir
    repoPaths, projectPaths = drop.setupPaths(projectRoot=projectDir)

    # create new directories
    for path in projectPaths.values():
        if not path.is_dir():
            path.mkdir(parents=True)
            logger.info(f"create {str(path)}")

    # hidden files
    copy_tree(str(wbuildPath), str(projectPaths["projectDir"] / ".wBuild"))
    copy_tree(str(repoPaths["modules"] / "helpers"), str(projectPaths["dropDir"] / "helpers"))
    # TODO: put version info there

    # copy Scripts and pipelines
    copyfile(repoPaths["template"] / "Snakefile", projectPaths["projectDir"] / "Snakefile")
    copy_tree(str(repoPaths["Scripts"]), str(projectPaths["Scripts"]))
    copyModuleCode(repoPaths, projectPaths)

    config_file = projectPaths["projectDir"] / "config.yaml"
    if not config_file.is_file():
        copyfile(repoPaths["template"] / "config.yaml", config_file)

    # search for a file containing the word readme and .md
    if not list(projectPaths["projectDir"].glob("readme*.md")):
        copyfile(repoPaths["template"] / "readme.md", projectPaths["projectDir"] / "readme.md")


@main.command()
def init():
    if Path(".drop").is_dir():
        print(".drop already exists, use drop update instead to update to a newer version")
    else:
        setFiles()
        logger.info("init...done")


@main.command()
def update():
    drop.checkDropVersion(Path().cwd().resolve(), force=True)
    logger.info("update...done")


@main.command()
def demo():
    # TODO: check if previous project gets overwritten
    setFiles()
    removeFile("config.yaml", warn=False)
    logger.info("init...done")

    # download data
    logger.info("download data")
    download_script = str(Path(drop.__file__).parent / "download_data.sh")
    response = subprocess.run(["bash", download_script], stderr=subprocess.STDOUT)
    response.check_returncode()

    # copy sample annotation and config files with absolute paths
    demo_repo = Path(drop.__file__).parent / "demo"
    drop.demo.fixSampleAnnotation(demo_repo / "sample_annotation_relative.tsv",
                                  Path.cwd() / "Data" / "sample_annotation.tsv")
    drop.demo.fixConfig(demo_repo / "config_relative.yaml", Path.cwd() / "config.yaml")

    logger.info("demo project created")
