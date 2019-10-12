import os
import oyaml

METHODS = {'AE': 'aberrant-expression-pipeline',
           'AS': 'aberrant-splicing-pipeline',
           'MAE': 'mae-pipeline'}
ROOT = os.path.join(os.getcwd(), ".drop")
#return os.path.join(os.path.dirname(__file__), "modules")

def getMethodPath(method, link_type='workdir', tmp_dir=None):
    """
    link_type: name of link flag for snakemake subworkflow
        workdir: directory of the submodule
        snakefile: path to Snakefile
        config_file: path to config file copy
        final_file: path to empty file used as last output of workflow
    """
    if method not in METHODS.keys():
        raise ValueError(f'{method} is not a valid method. Must be one of {METHODS.keys()}')
    
    if link_type == 'workdir':
        return os.path.join(ROOT, METHODS[method])
    elif link_type == 'snakefile':
        return os.path.join(ROOT, METHODS[method], 'Snakefile')
    elif link_type == 'rules':
        return os.path.join(ROOT, METHODS[method], 'snakeRules')
    elif link_type == 'config_file':
        assert tmp_dir is not None
        return os.path.join(tmp_dir, f'config_{method}.yaml')
    elif link_type == 'final_file':
        assert tmp_dir is not None
        return os.path.join(tmp_dir, f'{method}.done')
    else:
      raise ValueError(f'invalid link_type: "{link_type}"')

def setupTempFiles(config):
    
    # create temporary directory
    TMP_DIR = os.path.join(config['root'], 'tmp')
    if not os.path.exists(TMP_DIR):
        print(f"create temporary files directory {TMP_DIR}")
        os.mkdir(TMP_DIR)

    config_files = {}
    done_files = {}
    for method in METHODS.keys():
        
        # save config files
        conf_file = getMethodPath(method, link_type = 'config_file', tmp_dir=TMP_DIR)
        with open(conf_file, 'w') as yaml_file:
            oyaml.dump(config.copy(), yaml_file, default_flow_style=False)
            config_files[method] = conf_file
        
        # final rule output file
        done_file = getMethodPath(method, link_type='final_file', tmp_dir=TMP_DIR)
        done_files[method] = done_file
        # remove if it exists
        if os.path.exists(done_file):
            os.remove(done_file)

    return TMP_DIR, config_files, done_files

