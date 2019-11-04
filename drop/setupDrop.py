import drop
import subprocess
import pathlib
import wbuild

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
    print("install missing R packages")
    script = pathlib.Path(drop.__file__).parent / "installRPackages.R"
    requirements = str(pathlib.Path(drop.__file__).parent / 'requirementsR.txt')
    call = subprocess.Popen(
        ["Rscript", script, requirements], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT
    )
    
    stdout, stderr = call.communicate()
    if stderr:
        print(stderr)
        exit(1)
    stdout = stdout.decode()
    print(stdout)
    if "Execution halted" in stdout or "Error" in stdout:
        exit(1)

