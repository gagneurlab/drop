# -*- coding: utf-8 -*-
"""CLI interface to wbuild."""

import sys
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
    wbuildPath = pathlib.Path(wbuild.__file__).parent / '.wBuild'
    return templatePath, wbuildPath

@click.group()
@click_log.simple_verbosity_option(logger)
# @click.version_option('1.4.2',prog_name='wBuild')
def main():
    pass


@main.command()
def init():
    """Initialize the repository with wbuild.

    This will prepare wBuild in the current project
    """
    templatePath, wbuildPath = setup_paths()
    distutils.dir_util.copy_tree(str(wbuildPath), './.wBuild')
    if not os.path.isfile("Snakefile"):
        shutil.copy(str(templatePath / 'Snakefile'), '.')
    if not os.path.isfile("config.yaml"):
        shutil.copy(str(templatePath / 'config.yaml'), '.')
    if not os.path.exists("Scripts"):
        distutils.dir_util.copy_tree(str(templatePath), 'Scripts')
    
    ### search for a file containing the word readme and .md
    readme_exists = False
    onlyfiles = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
    for f in onlyfiles:
        if ("readme" in f) and f.endswith(".md"):
            readme_exists = True
            break
    if not readme_exists:
        #shutil.copy(str(templatePath / 'readme.md'), '.')
        open('readme.md', 'a').close()

    logger.info("init...done")