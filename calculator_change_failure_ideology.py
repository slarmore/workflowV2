#interface between the mol object and the actual calculators

import os
from shutil import Error
import subprocess
from . import molecule 
import numpy as np
from .utils import cleaner
from .message import warning,log,display


class Calculator:
    '''Calculator object is the middle-person between the actual 
    program and running the calculations'''
    
    def __init__(self,
                jobname,
                input,
                command,
                nproc,
                mem,
                time,
                partition,
                program,
                mol,
                argument_dict,
                try_count):

        #slurm inputs
        self.nproc = nproc
        self.mem = mem
        self.time = time
        self.partition = partition
        self.program = program
        self.jobname = jobname
        
        #file io
        self.dir = os.path.abspath(jobname) + '/'
        self.inputfile_full = self.dir +  self.jobname + 'try{0}'.format(try_count) + '.' + self.program.infiles[0]
        self.outputfile_full = self.dir +  self.jobname + 'try{0}'.format(try_count) + '.' + self.program.outfiles[0]
        self.inputfile_relative = self.jobname + 'try{0}'.format(try_count) + '.' + self.program.infiles[0]
        self.outputfile_relative = self.jobname + 'try{0}'.format(try_count) + '.' + self.program.outfiles[0]

        #input file
        self.input = input
        #write the input file in the right directory
        #creating the directory if it does not exist
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)

        for index,content in enumerate(self.input):
            with open('{0}/{1}-try{2}.{3}'.format(self.dir,self.jobname,self.try_count,self.program.infiles[index]),'w') as inputfile:
                inputfile.write(content)

        #submission 
        self.command = command
        #replace the generic command with the input/output files
        self.command = self.command.replace('INPUTFILE',self.inputfile_relative)
        self.command = self.command.replace('OUTPUTFILE',self.outputfile_relative)

        #always return a modified copy of the mol, rather than the modifying the mol itself
        self.mol = mol.copy()

        #resubmission
        self.try_count = try_count
        self.argument_dict = argument_dict

        #can be called by the Run or RunBatch functions to update the calculator
        #object based on the failures read from the current submission
        def resubmit(self):
            updated_argument_dict = self.program.fix_errors(self.mol,argument_dict)
            updated_calculator = self.program(updated_argument_dict,try_count=self.try_count+1)
            return(updated_calculator)

        #miscilaneous helper functions
        def InputFile(self):
            display(self.input)
        def Command(self):
            display(self.command)



def Run(calculator,tries=1):
    '''just submit a single job'''
    #check that the input is a calculator object
    if not isinstance(calculator,Calculator):
        raise TypeError('Expected Calcluator object')
    
    calculator.mol.warnings = ['First_submission']

    #submit as many times as is requested until the job finishes without warnings
    while len(calculator.mol.warnings) > 0:

        #if it's not the first submission, update the calculator 
        #with the resubmisson method
        if not calculator.mol.warnings[0] == 'First_submission':
            calculator = calculator.resubmit()

        if calculator.try_count <= tries:
            sbatch_name = '{0}-try{1}.sbatch'.format(calculator.jobname,calculator.try_count)
            with open('{0}{1}'.format(calculator.dir,sbatch_name),'w') as sbatch: #directory is already ended with /
                sbatch.write(generic_single(calculator))

            #submit the job
            log('submitting - sbatch {0}'.format(sbatch_name))
            p = subprocess.Popen('sbatch {0}'.format(sbatch_name), stdout=subprocess.PIPE, shell=True,cwd=calculator.dir)
            output, err = p.communicate()
            p_status = p.wait()
            slurmID = str(output).strip().split()[-1][:-3]
            slurmoutput = '{0}{1}-{2}.out'.format(calculator.dir,calculator.jobname,slurmID) #directory is already ended with /

            calculator.mol = calculator.program.read_output(calculator,slurmoutput)

        else:
            raise Error('Ran out of resubmission tries')

    return(calculator.mol)



def RunBatch(calculators,jobname='batch_job',max=50):
    '''Take a list of calculators and run all of them as a slurm array - assume they all have the same resources as the first one'''
    #check that the input is a calculator object
    for calculator in calculators:
        if not isinstance(calculator,Calculator):
            raise TypeError('Expected list of Calcluator objects')
    
    #create the directory to keep all of the batch job scripts
    if not os.path.isdir(jobname):
        os.makedirs(jobname)

    output = [None for calculator in calculators]

    #these are the indicies of the calculators that need to be submitted
    need_resub = list(range(0,len(calculators)))

    #while there are jobs with warnings, keep resubmitting
    while len(need_resub) > 0:

        #unless the number of tries is up
        if calculator.try_count <= calculator.tries:

            with open(jobname + '/' + jobname+ '.sbatch','w') as sbatch:
                sbatch.write(generic_batch(jobname,calculators,max,need_resub))

            for index in enumerate(need_resub):
                command_file_name = '{0}/{0}-{1}-try{2}.sh'.format(jobname,index+1,calculators[index].try_count)   #slurm array index start at 1 not 0
                with open(command_file_name,'w') as command_file:
                    command_file.write(generic_command(calculator))
        
            #submit the job
            log('submitting - sbatch {0}'.format(jobname))
            p = subprocess.Popen('sbatch {0}.sbatch'.format(jobname), stdout=subprocess.PIPE, shell=True,cwd=jobname)
            output, err = p.communicate()
            p_status = p.wait()
            slurmID = str(output).strip().split()[-1][:-3]

        
            #gather the output and account for failures
            new_need_resub = []
            for index in enumerate(need_resub):
                slurmoutput = '{0}/{0}-{1}_{2}.out'.format(jobname,slurmID,index+1) 
                calculators[index].mol = calculator[index].program.read_output(calculators[index],slurmoutput)
                output[index] = calculators[index].mol

                #if the returned mol object has warnings, 
                #use the calculator.resubmit() method to update the calculator
                if len(calculators[index].mol.warnings) > 0:
                    new_need_resub.append(index)
                    calculators[index].resubmit()

            need_resub = new_need_resub

        else:
            raise Error('Ran out of resubmission tries')
    
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

def generic_batch(jobname,calculators,max,need_resub):
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

""".format(jobname,partition,mem,nproc,time,len(need_resub),max)
    return(generic_batch_string)

def generic_command(calculator):
    generic_command_string = r"""#/bin/bash
#move to the directory with the input files
work={0}
cd $work

{1}
""".format(calculator.dir,calculator.command)
    return(generic_command_string)



