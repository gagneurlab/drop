from .setupDrop import setupDrop as drop
from .configHelper import ConfigHelper as config
from .submodules import *

def init():
    wbuild.cli.init()
    # compy our template

def update():
    wbuild.cli.update()

if __name__ == '__main__':
    import sys
    import wbuild
    
    arg = sys.args[1]
    if arg == 'init':
        init()
    elif arg == 'update':
        update()

