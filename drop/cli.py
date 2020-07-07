import wbuild
import drop
import yaml
from pathlib import Path
from shutil import copy
from distutils.dir_util import copy_tree, remove_tree
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
    return templatePath, modulePath, wbuildPath

def setFiles():
    templatePath, modulePath, wbuildPath = setup_paths()
    copy_tree(str(wbuildPath), "./.wBuild")
    
    copy(str(templatePath / "Snakefile"), ".")
    
    ### search for a file containing the word readme and .md
    if not len(list(Path(".").glob("readme*.md"))) > 0:
        open("readme.md", "a").close()
    
    # copy config file
    config_file = Path("config.yaml").resolve()
    if not config_file.is_file():
        copy(str(templatePath / 'config.yaml'), '.')
    
    # copy analysis scripts
    if not Path("Scripts").is_dir():
        copy_tree(str(templatePath / 'Scripts'), 'Scripts')
    
    # copy code for submodules
    drop_dir = Path(".drop").resolve()
    if drop_dir.is_dir():
        remove_tree(drop_dir)
        print("overwriting module scripts")
    copy_tree(str(modulePath), str(drop_dir / "modules"))
    
    # create symlinks for submodules
    drop_dir_link = drop_dir / "modules" / ".drop"
    drop_dir_link.symlink_to(drop_dir, target_is_directory=True)
    config_link = drop_dir / "config.yaml"
    config_link.symlink_to(config_file) 


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
    
    setFiles()
    
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
