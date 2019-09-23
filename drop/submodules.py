import os
import oyaml

def getRoot():
    return os.path.join(os.path.dirname(__file__), "modules")

methods = {'AE': 'aberrant-expression-pipeline',
           'AS': 'aberrant-splicing-pipeline',
           'MAE': 'mae-pipeline'}

def getMethodPath(method, link_type='workdir'):
    """
    link_type: name of link flag for snakemake subworkflow
        workdir: directory of the submodule
        snakeflie: path to Snakefile
        configfile: link to config file
        rules: link to extra rules
    """
    if method not in methods.keys():
        raise ValueError(f'{method} is not a valid method. Must be one of {methods.keys()}')
    
    if link_type == 'workdir':
        return os.path.join(getRoot(), methods[method])
    elif link_type == 'snakefile':
        return os.path.join(getRoot(), methods[method], 'Snakefile')
    elif link_type == 'rules':
      return os.path.join(getRoot(), methods[method], 'snakeRules')
    else:
      raise ValueError(f'invalid link_type: "{link_type}"')

def setupTempFiles(config):
    
    # create temporary directory
    tmpdir = os.path.join(config["ROOT"], config["DATASET_NAME"], 'tmp')
    if not os.path.exists(tmpdir):
      print(f"create temporary files directory {tmpdir}")
      os.mkdir(tmpdir)

    config_files = {}
    dummy_files = {}
    for method in methods.keys():
        # save config files
        conf_file = f'{tmpdir}/config_{method}.yaml'
        with open(conf_file, 'w') as yaml_file:
             oyaml.dump(config.copy(), yaml_file, default_flow_style=False)
             config_files[method] = conf_file
        # remove dummy files
        done_file = f'{tmpdir}/{method}.done'
        dummy_files[method] = done_file
        if os.path.exists(done_file):
            os.remove(done_file)

    return tmpdir, config_files, dummy_files

