#interface with X

########################################################################################################
#imports
from shutil import Error
from .. import calculator
from .. import molecule
from ..config import * #this is where the exe should be defined
import re
import numpy as np
from ..utils import cleaner
from ..message import warning,log,display
import os



########################################################################################################
#calculator creation function
def x(mol,jobname,runtype,method,nproc=1,mem=1,time=default_time,partition=default_partition,try_count=0,**kwargs):
    '''Create calculator object for X calculation'''
    

#######################
#check for valid input#
    if not isinstance(mol,molecule.Mol):
        if isinstance(mol,molecule.Conformer):
            #upgrade it to a mol object 
            mol = mol.ToMol()
        else:
            raise TypeError('Expected Mol or Conformer object')

    #always return a modified copy of the mol, rather than the modifying the mol itself
    mol = mol.copy()
    
    #key = runtype, value = list of acceptable arguments passed for calculation
    supported_runtypes = {

    }

    if not runtype in supported_runtypes:
        raise NotImplementedError('The {0} runtype is not implemented\n\nPlease use:\n{1}'.format(runtype,','.join([key for key in supported_runtypes])))


############################
#parse calculation keywords#

    #this is likely software spceific... but will likely need to loop over the passed keys
    for key, value in kwargs.items():
        if key in supported_runtypes[runtype]:

        #do something

        else:
            raise SyntaxError('The {0} keyword does not match the available options for the {1} runtype:\n{2}'.format(key,runtype,','.join(supported_runtypes[runtype])))


#########################   
#construct input file(s)#
    directory = os.path.abspath(jobname) + '/'
    basename_full = directory +  jobname + '-try{0}'.format(try_count) 
    
    #create string with the input file


########################
#define the run command#
    command = 'X INPUTFILE > OUTPUTFILE'.format(dir,program_exe)


##########################################################
#create calculator object to actually run the calculation#
    #need to combine all of the arguments into a single dict for failure resubmissions
    argument_dict = kwargs
    argument_dict['runtype'] = runtype
    argument_dict['time'] = time
    argument_dict['partition'] = partition
    argument_dict['nproc'] = nproc
    argument_dict['mem'] = mem
    argument_dict['jobname'] = jobname
    argument_dict['mol'] = mol

    return(calculator.Calculator(jobname=jobname,
                input=input,
                command=command,
                nproc=nproc,
                mem=mem,
                time=time,
                partition=partition,
                program=x(),
                mol=mol,
                argument_dict=argument_dict,
                try_count=try_count))



########################################################################################################
#software class
class X:

###################
#define attributes#
    def __init__(self):
        self.program_name = 'X'
        self.infiles = ['in']
        self.outfiles = ['out']
        self.normal_termination_line = -1   #where to look to see if calculation was successful
        self.normal_termination_string = 'Normal termination of X'   #what to look for
        self.unessesary_files = ['Gau*']

    #######################################
    #keywords to search for in output file#
        #key = string match in output file, value = function to read that section of the output
        self.keywords = {
    
        }

        #precompile all the regrex for efficiency's sake
        self.keys = [key for key in self.keywords]
        self.keys_regrexs = [re.compile(key) for key in self.keys]

####################
#read_output method#
    def read_output(self,calculator,slurmoutput):
        output = open(calculator.outputfile_full,'r')
        output_lines = output.read().splitlines()
        
        #clear the property dict and warning list
        calculator.mol.properties = {}
        calculator.mol.warnings = []
        #remove the unessesary leftover files
        #add the directory
        matches = [calculator.dir + match for match in self.unessesary_files]
        cleaner(matches)

        #for each line, check if any keywords are found
        #if so, use the matching function to update the mol object
        for line_number,line in enumerate(output_lines):
            for key,key_regrex in zip(self.keys,self.keys_regrexs):
                if key_regrex.search(line):
                    self.keywords[key](calculator.mol,line_number,line,output_lines,calculator)
        

        if not re.search(self.normal_termination_string,output_lines[self.normal_termination_line]):
            with open(slurmoutput,'r') as slurm:
                slurmoutput_content = slurm.read()

                warning('error in {0} calculation\nSlurm output:\n\n{1}\n\n'.format(self.program_name,slurmoutput_content))
            
            warning('last 10 lines of the output:\n\n{0}'.format('\n'.join(output_lines[-10:-1])))

            calculator.mol.warnings.append('{0}_Abnormal_Termination'.format(self.program_name))

            output.close()

        return(calculator.mol)

###################
#fix_errors method#
    def fix_errors(self,mol,input_name,kwargs):
        #do something
        return(kwargs)

#################
#resubmit method#
    def resubmit(self,kwargs,try_count):
        return(X(**kwargs,try_count=try_count))



########################################################################################################
#Functions for parsing output
#each of these functions should always take the same arguments
def get_some_value(mol,line_number,line,output_lines,calculator):
 
