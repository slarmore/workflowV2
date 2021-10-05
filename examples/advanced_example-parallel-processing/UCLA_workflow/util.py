import time
import pandas as pd
import random
import os
import shutil
import shutil
from workflowV2 import message
import multiprocessing as mp
import pickle
import queue
import time
import subprocess
import re

SENTINEL = None

def do_work(tasks_pending, tasks_completed):
    # Get the current worker's name
    worker_name = mp.current_process().name

    while True:
        try:
            task = tasks_pending.get_nowait()
        except queue.Empty:
            #message.log(worker_name + ' found an empty queue. Sleeping for a while before checking again...')
            time.sleep(0.01)
        else:
            try:
                if task == SENTINEL:
                    message.log(worker_name + ' no more work left to be done. Exiting...')
                    break

                time_start = time.perf_counter()
                work_func = pickle.loads(task['func'])
                message.log(worker_name + f' ({work_func.__name__}) received some work... ')
                result = work_func(**task['task'])
                tasks_completed.put({work_func.__name__: result})
                time_end = time.perf_counter() - time_start
                message.log(worker_name + f' ({work_func.__name__})done in {round(time_end,5)} seconds')
            except Exception as e:
                message.log(worker_name + f'({work_func.__name__}) task failed. ' + str(e))
                tasks_completed.put({work_func.__name__: None})


def par_proc(job_list, num_cpus=1):

    # Get the number of cores
    message.log('* Parallel processing')
    message.log('* Running on {} cores'.format(num_cpus))

    # Set-up the queues for sending and receiving data to/from the workers
    tasks_pending = mp.Queue()
    tasks_completed = mp.Queue()

    # Gather processes and results here
    processes = []
    results = []

    # Count tasks
    num_tasks = 0

    # Add the tasks to the queue
    for job in job_list:
        for task in job['tasks']:
            expanded_job = {}
            num_tasks = num_tasks + 1
            expanded_job.update({'func': pickle.dumps(job['func'])})
            expanded_job.update({'task': task})
            tasks_pending.put(expanded_job)

    # Use as many workers as there are cores (usually chokes the system so better use less)
    num_workers = num_cpus

    # We need as many sentinels as there are worker processes so that ALL processes exit when there is no more
    # work left to be done.
    for c in range(num_workers):
        tasks_pending.put(SENTINEL)

    message.log('* Number of tasks: {}'.format(num_tasks))

    # Set-up and start the workers
    for c in range(num_workers):
        p = mp.Process(target=do_work, args=(tasks_pending, tasks_completed))
        p.name = 'worker' + str(c)
        processes.append(p)
        p.start()

    # Gather the results
    completed_tasks_counter = 0
    while completed_tasks_counter < num_tasks:
        results.append(tasks_completed.get())
        completed_tasks_counter = completed_tasks_counter + 1

    for p in processes:
        p.join()

    return results
#a function to check for chk to read in and restart a calculation from


def formchk(file):
    p = subprocess.Popen(f'/work/lopez/g16/formchk {file}',stdout=subprocess.PIPE,shell=True,encoding='utf8')
    output,err = p.communicate()
    p_status = p.wait()
    
    #check if sucessful
    if re.search('Error termination',output):
        return(False)
    else:
        return(True)
    

def look_for_restart(mol_names,title,try_count=3):

    oldchk = [None] * len(mol_names)
    geom = [None] * len(mol_names)
    guess = [None] * len(mol_names)

    for index,name in enumerate(mol_names):
        while try_count > -1:
            chk_name = f"{name}-{title}/{name}-{title}-try{try_count}.chk"

            #check that the chk is actually readable
            if os.path.exists(chk_name) and formchk(chk_name):
                oldchk[index] = f"{name}-{title}/previous.chk"
                geom[index] = 'check'
                guess[index] = 'read'
                shutil.copyfile(chk_name,oldchk[index])
                message.log(f"Reading in previous {chk_name} as {oldchk[index]}")
                break
            try_count -= 1


    return(oldchk,geom,guess)

##################################################################

#I highly doubt two processes will finish computing at the same time
#and both go to write to the output csv... but in case they do,
#write out a file to lock access

#this is probably super unecessary now, but A. good practice
#and B. might be important as I scale up to more molecule/complexity

def open_csv(output_energies,who):

    #check whether it's free, if not wait a bit
    count = 0
    while os.path.exists(output_energies.split('.')[0]+'.lock'):
        #if we are waiting forever, assume that the .lock file was not deleted
        #nothing should take that long to write
        if count > 100:
            break
        message.log(f'{who} is waiting to read csv')
        time.sleep(random.randint(4,10))
        count += 1
        
    #once its free, lock it from other processes
    lock_csv(output_energies,who)

    return(pd.read_csv(output_energies))

def lock_csv(output_energies,who):
    with open(output_energies.split('.')[0]+'.lock','w') as lock:
        lock.write(f'{who} is locking {output_energies}')
    message.log(f'csv locked by {who}')

def unlock_csv(output_energies,who):
    os.remove(output_energies.split('.')[0]+'.lock')
    message.log(f'{who} unlocked csv')


def radical_electrons(mol):
    nradicals = 0

    for atom in mol.GetAtoms():
        nradicals += atom.GetNumRadicalElectrons()

    return(nradicals)



