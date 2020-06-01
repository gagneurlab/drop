import click
import wbuild
import drop
import yaml
import pathlib
import shutil
import distutils.dir_util
import subprocess
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
    templatePath = pathlib.Path(drop.__file__).parent / "template"
    modulePath = pathlib.Path(drop.__file__).parent / "modules"
    wbuildPath = pathlib.Path(wbuild.__file__).parent / ".wBuild"
    return templatePath, modulePath, wbuildPath

def setFiles():
    templatePath, modulePath, wbuildPath = setup_paths()
    distutils.dir_util.copy_tree(str(wbuildPath), "./.wBuild")
    
    shutil.copy(str(templatePath / "Snakefile"), ".")
    
    ### search for a file containing the word readme and .md
    if not len(list(pathlib.Path(".").glob("readme*.md"))) > 0:
        open("readme.md", "a").close()
    if not pathlib.Path("config.yaml").is_file():
        shutil.copy(str(templatePath / 'config.yaml'), '.')
    if not pathlib.Path("Scripts").is_dir():
        distutils.dir_util.copy_tree(str(templatePath / 'Scripts'), 'Scripts')
    
    if pathlib.Path(".drop").is_dir():
        distutils.dir_util.remove_tree(".drop")
        print("overwriting module scripts")
    distutils.dir_util.copy_tree(str(modulePath), ".drop/modules")


@main.command()
def init():
    if pathlib.Path(".drop").is_dir():
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
    
    setFiles()
    
    # download data
    logger.info("download data")
    download_script = str(pathlib.Path(drop.__file__).parent / "download_data.sh")
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
        dict_[key] = str(pathlib.Path(dict_[key]).resolve())
    
    with open("config.yaml", "w") as f:
        yaml.safe_dump(config.copy(), f, default_flow_style=False,
                       sort_keys=False)
    
    logger.info("demo project created")
