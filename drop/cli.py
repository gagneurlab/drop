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
logger = logging.getLogger(__name__)
click_log.basic_config(logger)

@click.group()
@click_log.simple_verbosity_option(logger)
@click.version_option('0.9.0',prog_name='drop')
def main():
    pass

def setup_paths():
    templatePath = Path(drop.__file__).parent / "template"
    modulePath = Path(drop.__file__).parent / "modules"
    wbuildPath = Path(wbuild.__file__).parent / ".wBuild"
    projectPath = Path.cwd().resolve()
    return projectPath, templatePath, modulePath, wbuildPath

def copyModuleCode(projectPath, templatePath, modulePath):
    # copy code for each submodule
    
    repo_map = {
        "aberrant-expression-pipeline": "AberrantExpression",
        "aberrant-splicing-pipeline": "AberrantSplicingAnalysis",
        "mae-pipeline": "MAEAnalysis"
    }
    
    helpers_path = projectPath/".drop"/"helpers"
    copy_tree(str(modulePath/"helpers"), str(helpers_path))
    
    for repo, analysis_dir in repo_map.items():
        repo_path = modulePath / repo
        workdir = projectPath / "Scripts" / analysis_dir / "pipeline"
        if workdir.is_dir():
            remove_tree(workdir)
            print(f"overwriting pipeline scripts for {analysis_dir}")
        copy_tree(str(repo_path), str(workdir))

def removeFile(filePath, warn=True):
    filePath = Path(filePath)
    if filePath.is_file():
        if warn:
            input(f"The file {str(filePath)} will be overwritten. Press Enter to continue...")
        filePath.unlink()
    
def setFiles(warn=True):
    projectPath, templatePath, modulePath, wbuildPath = setup_paths()
    
    # hidden files
    copy_tree(str(wbuildPath), ".wBuild")
    mkpath(".drop")
    # TODO: put version info there
    
    # copy Scripts and pipelines
    copyfile(templatePath/"Snakefile", projectPath/"Snakefile")
    copy_tree(str(templatePath/"Scripts"), str(projectPath/"Scripts"))
    copyModuleCode(projectPath, templatePath, modulePath)
    
    # config file
    config_file = Path("config.yaml")
    if not config_file.is_file():
        copyfile(templatePath/"config.yaml", config_file)
    
    # search for a file containing the word readme and .md
    if not len(list(projectPath.glob("readme*.md"))) > 0:
        open("readme.md", "a").close()


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
                 "genome": ["mae"], "qcVcf": ["mae"]}
    
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
