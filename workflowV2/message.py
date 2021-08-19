import sys
import os
from . import __logging__,__logfile__
from datetime import datetime

global __logging__
global __logfile__

def nologging(boolean):
    global __logging__
    if boolean:
        __logging__ = False
    else:
        __logging__ = True

def logtofile(file):
    global __logfile__
    global __logging__
 
    __logging__ = True

    if file is None:
        __logfile__ = None
    else:
        __logfile__ = os.path.abspath(file)

def warning(string,time=True):
    global __logging__
    global __logfile__

    if __logging__:

        if time:
            prefix = 'WORKFLOW WARNING:\n{0}\n    '.format(datetime.now())
        else:
            prefix = 'WORKFLOW WARNING:\n   '
        appendix = '\n\n'
        string = prefix + string + appendix

        if __logfile__ is None:
            sys.stderr.write(string)
            sys.stdout.flush()
        else:
            with open(__logfile__,'a') as log:
                log.write(string)

def log(string,time=True):
    global __logging__
    global __logfile__

    if __logging__:

        prefix = ''
        if time:
            appendix = ' - {0}\n'.format(datetime.now())
        else:
            appendix = '\n'
        string = prefix + string + appendix

        if __logfile__ is None:
            sys.stdout.write(string)
            sys.stdout.flush()
        else:
            with open(__logfile__,'a') as log:
                log.write(string)
        
def display(string):
    print(string,flush=True)
