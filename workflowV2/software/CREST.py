#interface with CREST

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
def CREST(mol,jobname,runtype,nproc=1,mem=1,time=default_time,partition=default_partition,max_confs=500,try_count=0,delete=['METADYN*','MRMSD','NORMMD*','*.tmp','wbo'],**kwargs):
    '''Create calculator object for CREST conformer generation'''
    

#######################
#check for valid input#
    if not isinstance(mol,molecule.Mol):
        if isinstance(mol,molecule.Conformer):
            #upgrade it to a mol object 
            mol = mol.ToMol()
        else:
            raise TypeError('Expected Mol or Conformer object')
    
    #key = runtype, value = list of acceptable arguments passed for calculation
    supported_runtypes = {'confsearch':['ewin','rthr','ethr','bthr','temp','esort','chrg','uhf',
                      'gbsa','alpb','opt','mdlen','shake','tstep','vbdump','zs','nocross',
                      'mrest','nci','tnmd','cbonds','cheavy','clight']
    }

    if not runtype in supported_runtypes:
        raise NotImplementedError('The {0} runtype is not implemented\n\nPlease use:\n{1}'.format(runtype,','.join([key for key in supported_runtypes])))


############################
#parse calculation keywords#

    #make sure the charge and multiplicity are correct
    if not 'chrg' in kwargs and mol.charge != 0:
        kwargs['chrg'] = mol.charge
    if not 'uhf' in kwargs and mol.mult != 1:
        kwargs['uhf'] = (mol.mult - 1 ) / 2

    arguments = ['-verbose']            #start with keywords always used
    for key, value in kwargs.items():
        if key in supported_runtypes[runtype]:
            if isinstance(value,bool):
                if value:                                  #if the option is set to True, include it, otherwise don't
                    arguments.append('-{0}'.format(key))
            else:                                                  #if it's not a bool, then assume it's setting a value
                arguments.append('-{0} {1}'.format(key,value))
        else:
            raise SyntaxError('The {0} keyword does not match the available options for the {1} runtype:\n{2}'.format(key,runtype,','.join(supported_runtypes[runtype])))

#########################   
#construct input file(s)#
    directory = os.path.abspath(jobname) + '/'
    basename_full = directory +  jobname + '-try{0}'.format(try_count) 

    #the input for CREST is just an xyz file, unless constraints are present, 
    # then need additional constraint file
    xyz = '{0}\n{1}\n{2}'.format(mol.natoms,jobname,mol.xyzstring)
    
    #create string with the input file
    #constraints if necessary
    if len(mol.constraints) > 0:
        if not '-subrmsd' in arguments:
            arguments.append('-subrmsd')
        arguments.append('-cinp constrain.c'.format(jobname,try_count))

        constraintfile = ['$constrain']
        
        constrained_atoms = []
        for constraint in mol.constraints:
            if len(constraint) == 2:
                constraintfile.append('distance: {0}, {1}, auto'.format(constraint[0]+1,constraint[1]+1))
                constrained_atoms.append(constraint[0]+1)
                constrained_atoms.append(constraint[1]+1)
            elif len(constraint) == 3:
                constraintfile.append('angle: {0}, {1}, {2}, auto'.format(constraint[0]+1,constraint[1]+1,constraint[2]+1))
                constrained_atoms.append(constraint[0]+1)
                constrained_atoms.append(constraint[1]+1)
                constrained_atoms.append(constraint[2]+1)
            elif len(constraint) == 4:
                constraintfile.append('dihedral: {0}, {1}, {2}, {3}, auto'.format(constraint[0]+1,constraint[1]+1,constraint[2]+1,constraint[3]+1))
                constrained_atoms.append(constraint[0]+1)
                constrained_atoms.append(constraint[1]+1)
                constrained_atoms.append(constraint[2]+1)
                constrained_atoms.append(constraint[3]+1)
        
        constraintfile.append('force constant={0}'.format(crest_constraint_force_constrant))
        constraintfile.append('reference=ref-try{0}.ref'.format(try_count))
        constraintfile.append('$metadyn')

        #get the list of atoms NOT constrained to include in the metadynamics
        missing = []
        constrained_atoms = list(set(constrained_atoms))
        constrained_atoms.sort()
        numbers = constrained_atoms
        numbers.insert(0,0) # add the minimum value on begining of the list
        numbers.append(mol.natoms+1)  # add the maximum value at the end of the list
        for rank in range(0, len(numbers)-1):
            if numbers[rank+1] - numbers[rank] > 2:
                missing.append("%s-%s"%(numbers[rank] +1 , numbers[rank+1] - 1))
            elif numbers[rank+1] - numbers[rank] == 2:
                missing.append(str(numbers[rank]+1))
        missing = str(missing)[1:-1]
        include = missing.replace("'","")

        constraintfile.append('atoms: {0}'.format(include))

        constraintfile.append('$end')
        constraintfile = '\n'.join(constraintfile)

    #there are a 3 CREST input files
    #1. the xyz
    #2. the constraint
    #3. the reference

    #the calculator will expect up to 3 files and unpack them with the extensions in order of:
    #.xyz
    #.c
    #.ref (note this is the same as the xyz, its just the reference for the constraints)
    if len(mol.constraints) > 0:
        input = [xyz,constraintfile,xyz]
        
    else:
        input = [xyz]



########################
#define the run command#
    
    #combine all of the parsed route line arguments
    arguments = ' '.join(arguments)

    #execution command

    command = 'ulimit -s unlimited\nexport OMP_STACKSIZE={0}G\nexport OMP_NUM_THREADS={1},1\n{2} INPUTFILE -xnam {3} -T {1} {4} > OUTPUTFILE'.format(mem,nproc,crest_exe,xtb_exe,arguments,jobname,try_count)


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
    argument_dict['max_confs'] = max_confs

    return(calculator.Calculator(jobname=jobname,
                input=input,
                command=command,
                nproc=nproc,
                mem=mem,
                time=time,
                partition=partition,
                program=crest(delete=delete),
                mol=mol,
                argument_dict=argument_dict,
                try_count=try_count))



########################################################################################################
#software class
class crest:

###################
#define attributes#
    def __init__(self,delete=['METADYN*','MRMSD','NORMMD*','*.tmp','wbo']):
        self.program_name = 'crest'
        self.infiles = ['{dir}{jobname}-try{try_count}.xyz','{dir}constrain.c','{dir}ref-try{try_count}.ref']
        self.outfiles = ['out']
        self.normal_termination_line = -1   #where to look to see if calculation was successful
        self.normal_termination_string = 'CREST terminated normally.'   #what to look for
        self.unessesary_files = delete

    #######################################
    #keywords to search for in output file#
        #key = string match in output file, value = function to read that section of the output
        self.keywords = {
                          'number of unique conformers for further calc': confs #wrap up getting all 
                                                                                #conformers here
                                                                                #crest is a bit of a special case
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
        raise NotImplementedError('No error handling implemented for CREST')
        return(kwargs)

#################
#resubmit method#
    def resubmit(self,kwargs,try_count):
        raise NotImplementedError('No resubmission implemented for CREST')
        return(CREST(**kwargs,try_count=try_count))



########################################################################################################
#Functions for parsing output
#each of these functions should always take the same arguments
def confs(mol,line_number,line,output_lines,calculator):
    #crest is a bit different from other programs, this is the only function to get all
    #of the conformer information and update the mol conformers in one
    energies = []
    conformers = []

    nconfs = int(line.split()[-1])
    max_conformers = calculator.argument_dict['max_confs']
    if nconfs > max_conformers:
                warning('CREST generated more than the max_conformers ({0}), only retaining the lowest {0}'.format(max_conformers))
                nconfs = max_conformers

    energies = [float(confline.split()[-1]) for confline in output_lines[line_number+1:line_number+1+nconfs]]
 
    with open('{0}/crest_conformers.xyz'.format(calculator.dir),'r') as conformerfile:
        conformerslines = conformerfile.read().splitlines()
        lines_per_xyz_block = mol.natoms + 2
        for conf_id in range(0,nconfs):
            start = (conf_id * lines_per_xyz_block) + 2   #don't care about the natom and title lines
            end = ( conf_id + 1 ) * lines_per_xyz_block
            rawxyz = conformerslines[start:end]

            coords = []
            atoms = []
            xyz = []

            for line in rawxyz:
                atom,x,y,z = line.split()
                x = float(x)
                y = float(y)
                z = float(z)
                if len(atom) > 1:
                    atom = atom[0] + atom[1:].lower()
                xyz.append([atom,x,y,z])
                atoms.append(atom)
                coords.append([x,y,z])
    
            xyz = np.array(xyz)
            atoms = np.array(atoms)
            coords = np.array(coords)

            conformers.append(molecule.Conformer(atoms=atoms,coords=coords,energy=energies[conf_id],tags={'origin':'CREST'},constraints=mol.constraints,charge=mol.charge,mult=mol.mult))
    
    #return the original mol object with the updated conformers
    mol.conformers = conformers
