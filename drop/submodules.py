import pathlib
import yaml
from snakemake.logging import logger

METHODS = {'AE': 'aberrant-expression-pipeline',
           'AS': 'aberrant-splicing-pipeline',
           'MAE': 'mae-pipeline'}
ROOT = pathlib.Path.cwd() / '.drop'
MODULES = ROOT / 'modules'
TMP_DIR = ROOT /'tmp'

def getTmpDir(str_=True):
    p = TMP_DIR
    if str_:
        p = str(p)
    return p

def getConfFile(str_=True):
    p = TMP_DIR / 'config.yaml'
    if str_:
        p = str(p)
    return p

def setupTempFiles(config):
    # create temporary directory
    if not TMP_DIR.exists():
        logger.info(f"create temporary files directory {TMP_DIR}")
        TMP_DIR.mkdir(parents=True)

    # save config file
    CONF_FILE = getConfFile(config)
    with open(CONF_FILE, 'w') as f:
        yaml.safe_dump(config.copy(), f)
    
    done_files = {}
    for method in METHODS.keys():
        
        # final rule output file
        done_file = getMethodPath(method, type_='final_file', str_=False)
        done_files[method] = str(done_file)
        
        # create module tmp Dir if missing
        tmp_dir = getMethodPath(method, type_='tmp_dir', str_=False)
        if not tmp_dir.exists():
            tmp_dir.mkdir(parents=True)
    
    return TMP_DIR, CONF_FILE, done_files

def getMethodPath(method, type_, str_=True):
    """
    type_: name of link flag for snakemake subworkflow
        workdir: directory of the submodule
        snakefile: path to Snakefile
        config_file: path to config file copy
        final_file: path to empty file used as last output of workflow
        unlock: path to empty file for unlocking subworkflow
    """
    if method not in METHODS.keys():
        raise ValueError(f'{method} is not a valid method. Must be one of {METHODS.keys()}')
    
    if type_ == 'workdir':
        p = MODULES / METHODS[method]
    elif type_ == 'snakefile':
        p = MODULES / METHODS[method] / 'Snakefile'
    elif type_ == 'tmp_dir':
        p = TMP_DIR / method
    elif type_ == 'final_file':
        p = TMP_DIR / f'{method}.done'
    elif type_ == 'unlock':
        p = TMP_DIR / f'{method}.unlock'
    else:
      raise ValueError(f'invalid type_: "{type_}"')
    
    if str_:
        p = str(p)
    return p


