import drop
import subprocess
import pathlib
import wbuild

def setupDrop(config):

    installRPackages()
    
    # parser config file
    parser = drop.config(config)
    config = parser.config
    config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)
    
    # temporary files
    tmpdir, config_files, dummy_files = drop.setupTempFiles(config)
    config["tmpdir"] = tmpdir
    config["configFileCopies"] = config_files
    config["finalFiles"] = dummy_files
    
    return parser, config

def callR(command):
    
    call = subprocess.Popen(
        ["R", "-e", f'"{command}"'], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT
    )
    
    stdout, stderr = call.communicate()
    if stderr:
        print(stderr)
    #print(stdout.decode())

def installRPackages():
    print("install missing R packages")
    command = """
    if (!requireNamespace('BiocManager', quietly = TRUE)) {
        install.packages('BiocManager')
    };

    packages <- readLines('requirementsR.txt');
    packages <- packages[packages != ''];
    for (package in packages) {
        if (!requireNamespace(package, quietly = TRUE):
        BiocManager::install(packages);
    };
    """
    
    callR(command)

