import wbuild
import drop
from pathlib import Path
from shutil import copy2, rmtree
from distutils.dir_util import copy_tree, remove_tree
import subprocess
import click
import click_log
import logging
import os
import filecmp as fc

wbuildPath = Path(wbuild.__file__).parent / ".wBuild"
logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.group()
@click_log.simple_verbosity_option(logger)
@click.version_option('1.3.1',prog_name='drop')


def main():
    pass


def overwrite(base_repo, local_proj):
    fc.clear_cache()  # clear file compare cache to avoid mistakes
    compare_obj = fc.dircmp(base_repo, local_proj)

    # remove all things not in the base_repo
    for i in compare_obj.right_only:
        logger.info(f"removing local file {i} it is not in the base drop")
        if os.path.isfile(local_proj / i):
            removeFile(local_proj / i, warn=False)
        else:
            remove_tree(local_proj / i)

    # for all dirs and files in base_dir
    for i in compare_obj.left_list:
        # files
        if os.path.isfile(base_repo / i):
            # filename is the same in both
            if i in compare_obj.common_files:

                # if file is diff copy original over. otherwise do nothing
                if i in compare_obj.diff_files:
                    logger.info(f"overwriting {local_proj / i} with {base_repo / i})")
                    copy2(base_repo / i, local_proj / i)


            # file not present in local project. Copy it
            else:
                logger.info(f"overwriting {local_proj / i} with {base_repo / i})")
                copy2(base_repo / i, local_proj / i)

        # dirs
        elif os.path.isdir(base_repo / i):
            if i in compare_obj.common_dirs:
                overwrite(base_repo / i, local_proj / i)
            else:
                logger.info(
                    f"the directory {str(base_repo / i)} does not exist locally. copying here: {str(local_proj)}")
                copy_tree(str(base_repo / i), str(local_proj / i))

        # other?
        else:
            logger.info(i, "is something other than file or dir. Ignoring")


def copyModuleCode(repoPaths, projectPaths):
    repo_map = {
        "aberrant-expression-pipeline": "AberrantExpression",
        "aberrant-splicing-pipeline": "AberrantSplicing",
        "mae-pipeline": "MonoallelicExpression",
        "rvc-pipeline": "rnaVariantCalling"

    }

    import sys
    for repo, analysis_dir in repo_map.items():
        fc.clear_cache()  # clear file compare cache to avoid mistakes

        base_repo = repoPaths["modules"] / repo
        local_proj = projectPaths["Scripts"] / analysis_dir

        #look for analysis_dir hidden from wbuild with "_" prefix and remove dir
        wbuild_hidden_path = projectPaths["Scripts"] / ("_" + analysis_dir)

        #if both hidden and local exist. Delete the hidden
        if wbuild_hidden_path.is_dir() and local_proj.is_dir():
            logger.info(f"removing the hidden wBuild path: {analysis_dir}")
            rmtree(wbuild_hidden_path,ignore_errors=True)
            logger.info("done")
        # if only hidden exists. rename and run normally
        elif wbuild_hidden_path.is_dir() and not local_proj.is_dir():
            logger.info(f"renaming the hidden wBuild path: {analysis_dir}")
            os.rename(wbuild_hidden_path,local_proj)
            logger.info("done")

        local_proj = local_proj / "pipeline"
        if not local_proj.is_dir():  # module directory does not exist. copy it
            logger.info(f"{local_proj} is not a directory, copy over from drop base")
            copy_tree(str(base_repo), str(local_proj))
        else:  # module dir does exist. Do a safe-overwrite
            logger.info(f"rewriting the module {analysis_dir} from the base DROP path")
            overwrite(base_repo, local_proj)
            logger.info("done")


def removeFile(filePath, warn=True):
    filePath = Path(filePath)
    if filePath.is_file():
        if warn:
            input(f"The file {str(filePath)} will be overwritten. Press Enter to continue...")
        filePath.unlink()


def setFiles(projectDir=None):
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
    copy2(repoPaths["template"] / "Snakefile", projectPaths["projectDir"] / "Snakefile")
    copy_tree(str(repoPaths["Scripts"]), str(projectPaths["Scripts"]))
    copyModuleCode(repoPaths, projectPaths)

    config_file = projectPaths["projectDir"] / "config.yaml"
    if not config_file.is_file():
        copy2(repoPaths["template"] / "config.yaml", config_file)

    # search for a file containing the word readme and .md
    if not list(projectPaths["projectDir"].glob("readme*.md")):
        copy2(repoPaths["template"] / "readme.md", projectPaths["projectDir"] / "readme.md")


@main.command()
def init():
    if Path(".drop").is_dir():
        print(".drop already exists, use drop update instead to update to a newer version")
    else:
        setFiles()
        logger.info("init...done")


@main.command()
def update():
    logger.info("updating local Scripts if necessary")
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
