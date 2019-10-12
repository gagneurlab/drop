import os
import click
import wbuild
import drop
import pathlib
import shutil
import distutils.dir_util
import click_log
import logging
logger = logging.getLogger(__name__)
click_log.basic_config(logger)

def setup_paths():
    """Setup the wbuild paths
    """
    templatePath = pathlib.Path(drop.__file__).parent / 'template'
    modulePath = pathlib.Path(drop.__file__).parent / 'modules'
    wbuildPath = pathlib.Path(wbuild.__file__).parent / '.wBuild'
    return templatePath, modulePath, wbuildPath

@click.group()
@click_log.simple_verbosity_option(logger)
@click.version_option('0.9.0',prog_name='drop')
def main():
    pass


def setFiles():
    templatePath, modulePath, wbuildPath = setup_paths()
    distutils.dir_util.copy_tree(str(wbuildPath), './.wBuild')
    
    if not os.path.isfile("Snakefile"):
        shutil.copy(str(templatePath / 'Snakefile'), '.')
    if not os.path.isfile("config.yaml"):
        shutil.copy(str(templatePath / 'config.yaml'), '.')
    if not os.path.exists("Scripts"):
        distutils.dir_util.copy_tree(str(templatePath), 'Scripts')
    
    if os.path.exists('.drop'):
        distutils.dir_util.remove_tree('.drop')
        print('overwriting module scripts')
    distutils.dir_util.copy_tree(str(modulePath), '.drop')


@main.command()
def init():
    
    if os.path.exists('.drop'):
        print('.drop already exists, use drop update instead to update to a newer version')
    else:
        setFiles()
        ### search for a file containing the word readme and .md
        readme_exists = False
        for f in os.listdir("."):
            if os.path.isfile(os.path.join(".", f)):
                continue
            if ("readme" in f) and f.endswith(".md"):
                readme_exists = True
                break
        if not readme_exists:
            open('readme.md', 'a').close()
        logger.info("init...done")

@main.command()
def update():
    # TODO: check version first
    setFiles()
    logger.info("update...done")
