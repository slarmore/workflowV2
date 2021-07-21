#interface with CREST

from . import calculator
from . import molecule
from .config import * #this is where the crest exe is defined
import re
import numpy as np
from .utils import cleaner
from .message import warning,log,display
from .config import max_conformers



def CREST(mol,jobname,runtype,nproc=1,mem=1,time='1-00:00:00',partition=default_partition,**kwargs):
    '''This function takes in a mol object and creates a calculator 
    object that can preform the requested calculation given 
    the setup here
    
    need to give 
    
    '''
    #check that the input is a Mol object
    if not isinstance(mol,molecule.Mol):
        raise TypeError('Expected Mol object')
    
    #always return a modified copy of the mol, rather than the modifying the mol itself
    mol = mol.copy()
    #supported CREST runtypes
    supported_runtypes = {

        'confsearch':['ewin','rthr','ethr','bthr','temp','esort','chrg','uhf',
                      'gbsa','alpb','opt','mdlen','shake','tstep','vbdump','zs','nocross',
                      'mrest','nci','tnmd']

                      #options that seem not to work no matter what I do... 'conds','cheavy','cmetal','clight'
    }
    
    if not runtype in supported_runtypes:
        raise NotImplementedError('The {0} runtype is not implemented for CREST\n\nPlease use:\n{1}'.format(runtype,','.join([key for key in supported_runtypes])))
        
    #parse all of these command line options
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

    #the input for CREST is just an xyz file, unless constraints are present, 
    # then need additional constraint file

    #xyz file
    xyz = '{0}\n{1}\n{2}'.format(mol.natoms,jobname,mol.xyzstring)

    #constraints if necessary
    if len(mol.constraints) > 0:
        if not '-subrmsd' in arguments:
            arguments.append('-subrmsd')
        arguments.append('-cinp {0}.c'.format(jobname))

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
        constraintfile.append('reference={0}.ref'.format(jobname))
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

    #combine all of the parsed route line arguments
    arguments = ' '.join(arguments)

    #execution command
    command = 'ulimit -s unlimited\nexport OMP_STACKSIZE={0}G\nexport OMP_NUM_THREADS={1},1\n{2} INPUTFILE -xnam {3} -T {1} {4} > OUTPUTFILE'.format(mem,nproc,crest_exe,xtb_exe,arguments)

    #return the calculator object to either run with calculator.Run or calculator.RunBatch
    return(calculator.Calculator(input=input,
                command=command,
                nproc=nproc,
                mem=mem,
                time=time,
                partition=partition,
                program=crest(),
                runtype=runtype,
                jobname=jobname,
                mol=mol))




class crest:
    def __init__(self):
        self.program_name = 'crest'
        self.infiles = ['xyz','c','ref']
        self.outfiles = ['out']

    def read_output(self,calculator,slurmoutput):
        #first check for normal termination
        output = open(calculator.outputfile_full,'r')
        output_lines = output.read().splitlines()
        
        #start at bottom of file where it's more likely to be found
        for line in reversed(output_lines):
            if re.search('CREST terminated normally.',line):

                #remove the unessesary leftover files
                matches = ['METADYN*','MRMSD','NORMMD*','*.tmp','wbo']
                #add the directory
                matches = [calculator.dir + match for match in matches]
                cleaner(matches)

                if calculator.runtype == 'confsearch':
                    result = confsearch(output_lines,calculator.dir,calculator.mol)
                output.close()
                return(result)
            else:
                with open(slurmoutput,'r') as slurm:
                    slurmoutput_content = slurm.read()

                    warning('error in CREST calculation\nSlurm output:\n\n{0}\n\n'.format(slurmoutput_content))
                
                warning('last 10 lines of the output:\n\n{0}'.format('\n'.join(output_lines[-10:-1])))

                output.close()
                return(None)

def confsearch(output_lines,dir,mol):
    '''Things to extract
            1. conformer energies
            2. conformer geometries - from the crest_conformers.xyz file


            needs to get the stuff required for a conformer object
            energy,atoms,coords,xyz=None,tags={},constraints=[]
    
    '''
    energies = []
    conformers = []


    for line_number,line in enumerate(output_lines):
        #line to look for vvv
        #number of unique conformers for further calc #number here 
        if re.search('number of unique conformers for further calc',line):
            nconfs = int(line.split()[-1])
            if nconfs > max_conformers:
                warning('CREST generated more than the max_conformers ({0}), only retaining the lowest {0}'.format(max_conformers))
                nconfs = max_conformers
                
            energies = [float(confline.split()[-1]) for confline in output_lines[line_number+1:line_number+1+nconfs]]
                

    with open('{0}/crest_conformers.xyz'.format(dir),'r') as conformerfile:
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
                if atom == 'CL':
                    atom = 'Cl'
                elif atom == 'BR':
                    atom = 'Br'
                xyz.append([atom,x,y,z])
                atoms.append(atom)
                coords.append([x,y,z])
    
            xyz = np.array(xyz)
            atoms = np.array(atoms)
            coords = np.array(coords)

            conformers.append(molecule.Conformer(energies[conf_id],atoms,coords,xyz=xyz,tags={'origin':'CREST'},constraints=mol.constraints))
    
    #return the original mol object with the updated conformers
    mol.conformers = conformers
    return(mol)




