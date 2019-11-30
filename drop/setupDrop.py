import drop
import subprocess
import pathlib
import wbuild
import re
from snakemake.logging import logger

def setupDrop(config):

    installRPackages()
    
    config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)
    parser = drop.config(config)

    tmp_dir, config_file, final_files = drop.setupTempFiles(parser.config)
    parser.setKey(parser.config, sub=None, key="tmpdir", default=tmp_dir)
    parser.setKey(parser.config, sub=None, key="configFile", default=config_file)
    parser.setKey(parser.config, sub=None, key="finalFiles", default=final_files)
    
    return parser, parser.parse()

def installRPackages():
    logger.info("check for missing R packages")
    script = pathlib.Path(drop.__file__).parent / "installRPackages.R"
    requirements = pathlib.Path(drop.__file__).parent / 'requirementsR.txt'
    
    #packages = [x.strip().split("#")[0] for x in open(requirements, 'r')]
    #packages = [x for x in packages if x != '']

    #for package in packages:
    #    logger.info(f"check {package}")   
    call = subprocess.Popen(
            ["Rscript", script, requirements], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT
        )
    
    # check output for errors
    stdout, stderr = call.communicate()
    if stderr:
        print(stderr)
        exit(1)
    stdout = stdout.decode()
    ep = re.compile("Execution halted|^ERROR", re.M)
    if ep.search(stdout):
        print(stdout)
        exit(1)

