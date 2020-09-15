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
@click.version_option('0.9.1', prog_name='drop')
def main():
    pass


def copyModuleCode(repoPaths, projectPaths):
    repo_map = {
        "aberrant-expression-pipeline": "AberrantExpression",
        "aberrant-splicing-pipeline": "AberrantSplicing",
        "mae-pipeline": "MonoallelicExpression"
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


def setFiles(warn=True):
    repoPaths, projectPaths = drop.setupPaths(projectRoot=Path.cwd().resolve())

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
    # TODO: check drop version first
    setFiles()
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

    # fix config file
    with open("config.yaml", "r") as f:
        config = yaml.load(f, Loader=yaml.Loader)
    path_keys = {"root": None,
                 "htmlOutputPath": None,
                 "sampleAnnotation": None,
                 "v29": ["geneAnnotation"],
                 "genome": None,
                 "qcVcf": ["mae"]}

    for key, sub in path_keys.items():
        # iterate to key and entry
        dict_ = config
        if sub is not None:
            for x in sub:
                dict_ = dict_[x]
        # set absolute path
        dict_[key] = str(Path(dict_[key]).resolve())

    with open("config.yaml", "w") as f:
        yaml.safe_dump(config.copy(), f, default_flow_style=False,
                       sort_keys=False)

    logger.info("demo project created")
