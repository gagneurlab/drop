import drop
import subprocess
import pathlib
import wbuild

def setupDrop(config):

    installRPackages()
    
    config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)
    parser = drop.config(config)

    tmpdir, confCopies, final = drop.setupTempFiles(parser.config)
    parser.setKey(parser.config, sub=None, key="tmpdir", default=tmpdir)
    parser.setKey(parser.config, sub=None, key="configFileCopies", default=confCopies)
    parser.setKey(parser.config, sub=None, key="finalFiles", default=final)
    
    return parser, parser.parse()

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

