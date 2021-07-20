#a collection of useful functions for the script

import os
import glob
import shutil
import multiprocessing as mp
from .message import warning,log,display
import datetime

atomic_number2symbol = {
    1:'H',
    6:'C',
    7:'N',
    8:'O',
    
}

symbol2atomic_number = {
    'H':1,
    'C':6,
    'N':7,
    'O':8,
    
}

#conversions
def hartree2kcal(hartree):
    return(hartree * 627.509)



def cleaner(matches):
    for match in matches:
        files = glob.glob(match)
        for file in files:
            if os.path.isdir(file):
                shutil.rmtree(file)
            else:
                os.remove(file)





def par_update(result):
    global results
    global complete
    global total

    i,item,result = result
    
    results[i] = result
    complete += 1
    log('Done processing {0} - {1}'.format(item,datetime.datetime.now()))
    log('{0:.2f}% complete'.format(complete/total))

def par_error(error):
    warning('Error in parallel processing: {0} - {1}'.format(error,datetime.datetime.now()))

def par_wrapper(i,func,item,**kwargs):
    log('Staring par processing {0} - {1}'.format(item,datetime.datetime.now()))
    return(i,item,func(item,**kwargs))

def parallelize(items,func,nproc,**kwargs):
    '''interface for parallelizing calculations on a list of mol objects'''
    global results
    global complete
    global total

    complete = 0
    total = len(items)
    results = [None] * total

    pool = mp.Pool(nproc)

    for i,item in enumerate(items):
        args = (i,func,item)
        pool.apply_async(par_wrapper,args,kwargs,callback=par_update, error_callback=par_error)

    pool.close()
    pool.join()
    return(results)




