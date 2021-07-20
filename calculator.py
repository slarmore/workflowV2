#interface between the mol object and the actual calculators

import os
import subprocess
from . import molecule 
import numpy as np
from .utils import cleaner
from .message import warning,log,display


class Calculator:
    ''' calculator object is created and serves as the IO interface
    
    #create the mol object
    mol = SmilesToMol('CCC')

    #create the calculator object by calling the setup for
    #whatever program is actually run
    gaussian_calc = Gaussian.opt(mol,params)

    #the the calculator object manages the writing files and submission
    result_mol = gaussian_calc.run()

    #that way I can do a batch submission 


    calculators = [Gauss.opt(conf) for conf in mol.conformers]
    RunParallel(calculators)
    
    to run them as an array job


    so need a 'Run','RunParallel', 'RunSerial' function to
    take in a calculator and write an appropriate sbatch

    
    so what does this object actually consist of?
    i think its mostly just a container, and then writes the files



    it should take in a string that is the input file,
    it should know the number of cores/memory/slurm stuff
    it should know the command to run
    if should know the program called

    '''

    
    def __init__(self,
                 input=None,
                 command=None,
                 nproc=None,
                 mem=None,
                 time=None,
                 partition=None,
                 program=None,
                 runtype=None,
                 jobname=None,
                 mol=None):


        self.input = input
        self.command = command
        self.nproc = nproc
        self.mem = mem
        self.time = time
        self.partition = partition
        self.program = program
        self.runtype = runtype
        self.jobname = jobname
        self.dir = os.path.abspath(jobname) + '/'
        self.inputfile_full = self.dir +  self.jobname + '.' + self.program.infiles[0]
        self.outputfile_full = self.dir +  self.jobname + '.' + self.program.outfiles[0]
        self.inputfile_relative = self.jobname + '.' + self.program.infiles[0]
        self.outputfile_relative = self.jobname + '.' + self.program.outfiles[0]
        self.mol = mol

        #write the input file in the right directory
        #creating the directory if it does not exist
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        
        for index,content in enumerate(self.input):
            with open('{0}/{1}.{2}'.format(self.dir,self.jobname,self.program.infiles[index]),'w') as inputfile:
                inputfile.write(content)

        #replace the generic command with the input/output files
        self.command = self.command.replace('INPUTFILE',self.inputfile_relative)
        self.command = self.command.replace('OUTPUTFILE',self.outputfile_relative)

        def InputFile(self):
            display(self.input)
        def Command(self):
            display(self.command)


def Run(calculator):
    '''just submit a single job'''
    #check that the input is a calculator object
    if not isinstance(calculator,Calculator):
        raise TypeError('Expected Calcluator object')

    sbatch_name = '{0}.sbatch'.format(calculator.jobname)
    with open('{0}{1}'.format(calculator.dir,sbatch_name),'w') as sbatch: #directory is already ended with /
        sbatch.write(generic_single(calculator))
    
    #submit the job
    log('submitting - sbatch {0}'.format(sbatch_name))
    p = subprocess.Popen('sbatch {0}'.format(sbatch_name), stdout=subprocess.PIPE, shell=True,cwd=calculator.dir)
    output, err = p.communicate()
    p_status = p.wait()
    slurmID = str(output).strip().split()[-1][:-3]
    slurmoutput = '{0}{1}-{2}.out'.format(calculator.dir,calculator.jobname,slurmID) #directory is already ended with /

    return(calculator.program.read_output(calculator,slurmoutput))
    

def RunBatch(calculators,jobname='batch_job',max=50):
    '''Take a list of calculators and run all of them as a slurm array - assume they all have the same resources as the first one'''
    #check that the input is a calculator object
    for calculator in calculators:
        if not isinstance(calculator,Calculator):
            raise TypeError('Expected list of Calcluator objects')
    
    #create the directory to keep all of the batch job scripts
    if not os.path.isdir(jobname):
        os.makedirs(jobname)
    
    with open(jobname + '/' + jobname+ '.sbatch','w') as sbatch:
        sbatch.write(generic_batch(jobname,calculators,max))
    for index,calculator in enumerate(calculators):
        command_file_name = '{0}/{0}-{1}.sh'.format(jobname,index+1)   #slurm array index start at 1 not 0
        with open(command_file_name,'w') as command_file:
            command_file.write(generic_command(calculator))
    
    #submit the job
    log('submitting - sbatch {0}'.format(jobname))
    p = subprocess.Popen('sbatch {0}.sbatch'.format(jobname), stdout=subprocess.PIPE, shell=True,cwd=jobname)
    output, err = p.communicate()
    p_status = p.wait()
    slurmID = str(output).strip().split()[-1][:-3]

    output = []
    for index,calculator in enumerate(calculators):
        slurmoutput = '{0}/{0}-{1}_{2}.out'.format(jobname,slurmID,index+1) #directory is already ended with /
        output.append(calculator.program.read_output(calculator,slurmoutput))
    return(output)


#############################################################################################################################
#templates
#############################################################################################################################

def generic_single(calculator):
    generic_single_string = r"""#!/bin/bash
#SBATCH --wait  #wait until the job finishes to release the shell - this will allow python to wait for the calculation end
#SBATCH --job-name={0}
#SBATCH --out={0}-%A.out
#SBATCH --partition={1}
#SBATCH --mem={2}000
#SBATCH --ntasks={3}
#SBATCH --nodes=1
#SBATCH --time={4}

#execute the calculation command
{5}

""".format(calculator.jobname,calculator.partition,calculator.mem,calculator.nproc,calculator.time,calculator.command)
    return(generic_single_string)

def generic_batch(jobname,calculators,max):
    #need to get the slurm parameters - get it from the first calculator and just assume they are the same....
    partition = calculators[0].partition
    mem = calculators[0].mem
    nproc = calculators[0].nproc
    time = calculators[0].time

    generic_batch_string = r"""#!/bin/bash
#SBATCH --wait  #wait until the job finishes to release the shell - this will allow python to wait for the calculation end
#SBATCH --job-name={0}
#SBATCH --out={0}-%A_%a.out
#SBATCH --partition={1}
#SBATCH --mem={2}000
#SBATCH --ntasks={3}
#SBATCH --nodes=1
#SBATCH --time={4}
#SBATCH --array=1-{5}%{6}

#keep each calculator command in a separate script that get's called
#based on the array task id
chmod 777 {0}-${{SLURM_ARRAY_TASK_ID}}.sh

./{0}-${{SLURM_ARRAY_TASK_ID}}.sh

""".format(jobname,partition,mem,nproc,time,len(calculators),max)
    return(generic_batch_string)

def generic_command(calculator):
    generic_command_string = r"""#/bin/bash
#move to the directory with the input files
work={0}
cd $work

{1}
""".format(calculator.dir,calculator.command)
    return(generic_command_string)



