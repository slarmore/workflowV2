#interface with GAUSSIAN
from shutil import Error
from .. import calculator
from .. import molecule
from ..config import * #this is where the g16 exe is defined
import re
import numpy as np
from ..utils import cleaner
from ..message import warning,log,display
import os


def GAUSSIAN(mol,jobname,runtype,method,nproc=1,mem=1,time='1-00:00:00',partition=default_partition,oldchk=None,**kwargs):
    '''This function takes in a mol object and creates a calculator 
    object that can preform the requested calculation given 
    the setup here
    
    need to give 
    
    '''
    #check that the input is a Mol or conformer type object
    if not isinstance(mol,molecule.Mol):
        if isinstance(mol,molecule.Conformer):
            #upgrade it to a mol object 
            mol = mol.ToMol()
        else:
            raise TypeError('Expected Mol or Conformer object')

    #always return a modified copy of the mol, rather than the modifying the mol itself
    mol = mol.copy()
    #supported Gaussian runtypes
    supported_runtypes = {
        'opt':['opt','IOp','method','scrf','scf','guess','pseudo','temperature','pop','density','afterinput','empiricaldispersion','geom','td'],
        'opt_freq':['opt','freq','IOp','method','scrf','scf','guess','pseudo','temperature','pop','density','afterinput','empiricaldispersion','geom','td'],
        'freq':['freq','IOp','method','scrf','scf','guess','pseudo','temperature','pop','density','afterinput','empiricaldispersion','geom','td'],
        'sp':['IOp','method','scrf','scf','guess','pseudo','temperature','pop','density','afterinput','empiricaldispersion','geom','td'],
        'irc':['irc','IOp','method','scrf','scf','guess','pseudo','temperature','pop','density','afterinput','empiricaldispersion','geom','td'],
        'tddft':['td','IOp','method','scrf','scf','guess','pseudo','temperature','pop','density','afterinput','empiricaldispersion','geom','opt']
    }

    if not runtype in supported_runtypes:
        raise NotImplementedError('The {0} runtype is not implemented for Gaussian\n\nPlease use:\n{1}'.format(runtype,','.join([key for key in supported_runtypes])))

    #add the modredundant section with the constraints 

    if len(mol.constraints) > 0:
        constraint_string = []
        for constraint in mol.constraints:
            if len(constraint) == 1:
                constraint_string.append('X {0} F'.format(constraint))
            elif len(constraint) == 2:
                constraint_string.append('B {0} {1} F'.format(constraint[0],constraint[1]))
            elif len(constraint) == 3:
                constraint_string.append('A {0} {1} {2} F'.format(constraint[0],constraint[1],constraint[2]))
            elif len(constraint) == 4:
                constraint_string.append('D {0} {1} {2} {3} F'.format(constraint[0],constraint[1],constraint[2],constraint[3])) 

        constraint_string = '\n'.join(constraint_string)

        if runtype in ['opt','opt_freq']:
            if 'opt' in kwargs:
                current = kwargs['opt']
                if not re.search('modredundant',current):
                    kwargs['opt'] = current + ' modredundant'
            else:
                kwargs['opt'] = 'modredundant'

        if 'afterinput' in kwargs:
            kwargs['afterinput'] += constraint_string
        else:
            kwargs['afterinput'] = constraint_string


    #if no arguments for the base runtype, set to default
    maintypes = runtype.split('_')               #expects multicomponent jobs to be called by STEP1_STEP2 (ie opt_freq)
    runcommands = []
    for maintype in maintypes:
        if not maintype in kwargs:
            if maintype == 'sp':               #Gaussian is horrible... sp() is not valid, but opt=() is.. so if doing a sp
                                               #need to just exclude the sp all togeher
                runcommands.append('')
            else:
                runcommands.append('{0}=()'.format(maintype))
        else:
            runcommands.append('{0}=({1})'.format(maintype,kwargs[maintype]))

    #parse all peripherial arguments - not the main runtype
    afterinputsection = []
    route = []
    for key, value in kwargs.items():
        if key in maintypes:
            pass
        elif key in supported_runtypes[runtype]:
            if key == 'afterinput':
                afterinputsection.extend([value,''])
            else:
                route.append('{0}=({1})'.format(key,value))
        else:
            raise SyntaxError('The {0} keyword does not match the available options for the {1} runtype:\n{2}'.format(key,runtype,','.join(supported_runtypes[runtype])))
    
    #combine all of the parsed route line arguments
    route = ' '.join(route)

    #to ensure that IOps will work, need to split up optfreq jobs into two steps
    if runtype == 'opt_freq':
        com = []
        #only include oldchk on the optimization step 
        if oldchk is not None:
            com.append('%oldchk={0}'.format(os.path.abspath(oldchk)))
        
        com.extend(['%chk={0}.chk'.format(jobname),
                '%nprocs={0}'.format(nproc),
                '%mem={0}GB'.format(mem),
                '# {0} {1} {2}'.format(method,runcommands[0],route),      #opt will be the first run command in the list
                '',
                jobname,
                '',
                '{0} {1}'.format(mol.charge,mol.mult)])

        if not 'geom' in kwargs:
            com.append(mol.xyzstring),
        com.append('')

        if len(afterinputsection) > 0:
            com.append('\n'.join(afterinputsection))
        
        #start the freq step input
        com.append('--Link1--')

        #add geom=check guess=read if not already present in the original route
        print('this is the route',route)
        if not re.search('geom=\(check\)',route):
            route += 'geom=(check) '
        if not re.search('guess=\(read\)',route):
            route += 'guess=(read) '

        com.extend(['%chk={0}.chk'.format(jobname),
                '%nprocs={0}'.format(nproc),
                '%mem={0}GB'.format(mem),
                '# {0} {1} {2}'.format(method,runcommands[1],route),    #freq will be the second run command in the list
                '',
                jobname,
                '',
                '{0} {1}'.format(mol.charge,mol.mult),
                ''])
        if len(afterinputsection) > 0:
            com.append('\n'.join(afterinputsection))
        com.append('')

    #if it is not opt_freq, then only need one input section    
    else:
        com = []
        #only include oldchk on the optimization step 
        if oldchk is not None:
            com.append('%oldchk={0}'.format(os.path.abspath(oldchk)))
        
        com.extend(['%chk={0}.chk'.format(jobname),
                '%nprocs={0}'.format(nproc),
                '%mem={0}GB'.format(mem),
                '# {0} {1} {2}'.format(method,' '.join(runcommands),route),    #concatenate the list of run commands
                '',
                jobname,
                '',
                '{0} {1}'.format(mol.charge,mol.mult)])

        if not 'geom' in kwargs:
            com.append(mol.xyzstring),
        com.append('')

        if len(afterinputsection) > 0:
            com.append('\n'.join(afterinputsection))
        com.append('')

    input = ['\n'.join(com)]

    #execution command
    command = 'export GAUSS_SCRDIR=./\nexport g16root={1}\n. $g16root/g16/bsd/g16.profile\n$g16root/g16/g16 INPUTFILE > OUTPUTFILE'.format(dir,g16root)

    #return the calculator object to either run with calculator.Run or calculator.RunBatch
    return(calculator.Calculator(input=input,
                command=command,
                nproc=nproc,
                mem=mem,
                time=time,
                partition=partition,
                program=gaussian(),
                runtype=runtype,
                jobname=jobname,
                mol=mol))


def GAUSSIAN_fixerrors(mol,jobname,runtype,method,nproc=1,mem=1,time='1-00:00:00',partition=default_partition,oldchk=None,maxtries=3,TS=False,**kwargs):
    mol.warnings = ['first_submission']
    tries = 0
    while len(mol.warnings) > 0:
        tries += 1

        #make fixes based on the warnings in mol.warnings
        if tries > 1:
            kwargs = fix_errors(mol,kwargs,TS)

        log('{0} submission {1}'.format(jobname,tries))
        calculator_obj = GAUSSIAN(mol,jobname='{0}-{1}'.format(jobname,tries),runtype=runtype,method=method,nproc=nproc,mem=mem,time=time,partition=partition,oldchk=oldchk,**kwargs)
        oldchk = '{0}-{1}/{0}-{1}.chk'.format(jobname,tries)

        if tries > maxtries:
            raise Error('Ran out of fixerror tries')

        mol = calculator.Run(calculator_obj)

def GAUSSIAN_fixerrors_batch(mols,jobname,runtype,method,nproc=1,mem=1,time='1-00:00:00',partition=default_partition,oldchk=None,maxtries=3,TS=False,**kwargs):
    if oldchk is None:
        oldchk = [None for mol in mols]
    
    else:
        if len(oldchk) != len(mols):
            raise IndexError('The same number of molecules and oldchk files must be given. Use None for mols that do not need oldchks')
    
    for mol in mols:
        mol.warnings = ['first_submission']

    tries = 0 

    need_submissions = [mol for mol in mols]

    while len(need_submissions) > 0:
        tries += 1

        #make fixes based on the warnings in mol.warnings
        if tries > 1:
            kwargs = [fix_errors(mol,kwargs[index],TS) if mol in need_submissions else kwargs[index] for index,mol in enumerate(mols)]
        else:
            kwargs = [kwargs.copy() for mol in mols]

        log('{0} submission {1}'.format(jobname,tries))
        calculators = [GAUSSIAN(mol,jobname='{0}-{1}-try{2}'.format(jobname,index,tries),runtype=runtype,method=method,nproc=nproc,mem=mem,time=time,partition=partition,oldchk=oldchk[index],**kwargs[index]) for index,mol in enumerate(mols)]
        oldchk = ['{0}-{1}-try{2}/{0}-{1}-try{2}.chk'.format(jobname,index,tries) for index,mol in enumerate(mols)]

        if tries > maxtries:
            raise Error('Ran out of fixerror tries')
        
        need_submissions = []
        for mol in mols:
            if len(mol.warnings) > 0:
                need_submissions.append(mol)

        mol = calculator.RunBatch(calculators)



def fix_errors(mol,kwargs,TS):

    #not very elegant, but loop through the list and 
    #add keywords to fix any known errors
    
    if 'Non-Optimized' in mol.warnings:
        log('Fixed Non-Optimized Error')
        if 'freq' in kwargs:
            if 'opt' in kwargs:
                current = kwargs['opt']
                #need a way to remove 'calcfc' but not remove 'recalcfc'
                if not re.search('readfc',current):
                    kwargs['opt'] = current + ',readfc'
            else:
                kwargs['opt'] = 'readfc'

            kwargs['geom'] = 'check'
            kwargs['guess'] = 'read'

        else:
            if 'opt' in kwargs:
                current = kwargs['opt']
                if not re.search('calcfc',current):
                    kwargs['opt'] = current + ',calcfc'
            else:
                kwargs['opt'] = 'calcfc'
    
    if 'SCF_Error' in mol.warnings:
        log('Fixed SCF_Error')
        if 'scf' in kwargs:
            current = kwargs['scf']
            if not re.search('qc'):
                kwargs['scf'] == current + ',qc'
        else:
            kwargs['scf'] = 'qc'
    
    if 'negative_frequency' in mol.warnings:
        if TS:
        #check for saddle point
            if 'frequencies' in mol.properties:
                if mol.properties['frequencies'][1] < 0:
                    log('Fixed saddle point')
                    if 'opt' in kwargs:
                        current = kwargs['opt']
                        if not re.search('readfc',current):
                            kwargs['opt'] = current + ',readfc'
                
                    kwargs['geom'] = 'check'
                    kwargs['guess'] = 'read'
    
        else:
        #check for TS
            log('Fixed negative frequency')
            if 'opt' in kwargs:
                current = kwargs['opt']
                if not re.search('readfc',current):
                    kwargs['opt'] = current + ',readfc'
            else:
                kwargs['opt'] = 'readfc'

            kwargs['geom'] = 'check'
            kwargs['guess'] = 'read'

    return(kwargs)



class gaussian:
    def __init__(self):
        self.program_name = 'gaussian'
        self.infiles = ['com']
        self.outfiles = ['log']

        #need a list of keywords to look for when parsing the output file
        #each keyword will point to a function that extracts the information from that keyword
        #and updates the mol object
        self.keywords = {
                             'SCF Done:': electronic_energies,
                             'Alpha  occ. eigenvalues': orbital_energies,
                             'Thermal correction to Energy': thermal_energies,
                             'Standard orientation': geometry,
                             'Non-Optimized Parameters': non_opt,
                             'armonic frequencies': frequencies,       #want to match both Harmonic and Anharmonic
                             'SCF Error SCF Error SCF Error SCF Error': scf_error,
                             'Total Energy, E': tddft_energy

        }

        self.keys = [key for key in self.keywords]
        self.keys_regrexs = [re.compile(key) for key in self.keys]


    def read_output(self,calculator,slurmoutput):
        global step 
        step = -1

        #first check for normal termination
        output = open(calculator.outputfile_full,'r')
        output_lines = output.read().splitlines()
        
        #clear the property dict and warning list
        calculator.mol.properties={}
        calculator.mol.warnings = []
        #remove the unessesary leftover files
        matches = ['Gau*']
        #add the directory
        matches = [calculator.dir + match for match in matches]
        cleaner(matches)
        #mol = calculator.mol

        #for each line, check if any keywords are found
        #if so, use the matching function to update the mol object
        for line_number,line in enumerate(output_lines):
            for key,key_regrex in zip(self.keys,self.keys_regrexs):
                if key_regrex.search(line):
                    self.keywords[key](calculator.mol,line_number,line,output_lines,calculator)
        

        if re.search('Normal termination of Gaussian',output_lines[-1]):
            return(calculator.mol)
        else:
            with open(slurmoutput,'r') as slurm:
                slurmoutput_content = slurm.read()

                warning('error in Gaussian calculation\nSlurm output:\n\n{0}\n\n'.format(slurmoutput_content))
            
            warning('last 10 lines of the output:\n\n{0}'.format('\n'.join(output_lines[-10:-1])))

            calculator.mol.warnings.append('Gaussian_Abnormal_Termination')

            output.close()
            return(calculator.mol)


def electronic_energies(mol,line_number,line,output_lines,calculator):
    energy = float(line.split()[4])
    mol.energy = energy
    
    #add to list of optimization energies for
    if 'optimization_energies' in mol.properties:
        mol.properties['optimization_energies'].append(energy)
    else:
        mol.properties['optimization_energies'] = [energy]


def orbital_energies(mol,line_number,line,output_lines,calculator):
    ''' put all of the occupied orbitals into the occupied_orbital_energies thing and all the unoccupied in its own'''

    #there will be multiple matches, so check that the next line contains the lumo, if not
    #then pass
    if re.search('Alpha virt. eigenvalues',output_lines[line_number+1]):
        occupied_orbitals = []
        homo_line_number = line_number
        homo_line = output_lines[homo_line_number].split()[4:]
        occupied_orbitals.extend(list(reversed([float(energy) for energy in homo_line])))

        #need to walk backward through the file and grab all of the occupied orbital energies 
        # until the line does not contain 'Alpha occ. eigenvalues'
        homo_line_number -= 1
        newline = output_lines[homo_line_number]
        while re.search('Alpha  occ. eigenvalues',newline):
            newline = newline.split()[4:]
            occupied_orbitals.extend(list(reversed([float(energy) for energy in newline])))
            homo_line_number -= 1
            newline = output_lines[homo_line_number]

        unoccupied_orbitals = []
        lumo_line_number = line_number + 1
        lumo_line = output_lines[lumo_line_number].split()[4:]
        unoccupied_orbitals.extend([float(energy) for energy in lumo_line])

        #need to walk forward through the file and grab all of the unoccupied orbital energies 
        # until the line does not contain 'Alpha virt. eigenvalues'
        lumo_line_number += 1
        newline = output_lines[lumo_line_number]
        while re.search('Alpha virt. eigenvalues',newline):
            newline = newline.split()[4:]
            unoccupied_orbitals.extend([float(energy) for energy in newline])
            lumo_line_number += 1
            newline = output_lines[homo_line_number]

        homo = occupied_orbitals[0]
        lumo = unoccupied_orbitals[0]

        mol.properties['homo'] = homo
        mol.properties['lumo'] = lumo
        mol.properties['occupied_orbitals'] = occupied_orbitals
        mol.properties['unoccupied_orbitals'] = unoccupied_orbitals
    
def thermal_energies(mol,line_number,line,output_lines,calculator):
    thermal_correction = float(line.split()[-1])
    enthalpy_correction = float(output_lines[line_number+1].split()[-1])
    free_energy_correction = float(output_lines[line_number+2].split()[-1])
    zero_point_corrected_energy = float(output_lines[line_number+3].split()[-1])
    thermal_energy = float(output_lines[line_number+4].split()[-1])
    enthalpy = float(output_lines[line_number+5].split()[-1])
    free_energy = float(output_lines[line_number+6].split()[-1])

    mol.properties['thermal_correction'] = thermal_correction
    mol.properties['enthalpy_correction'] = enthalpy_correction
    mol.properties['free_energy_correction'] = free_energy_correction
    mol.properties['zero_point_corrected_energy'] = zero_point_corrected_energy
    mol.properties['thermal_energy'] = thermal_energy
    mol.properties['enthalpy'] = enthalpy
    mol.properties['free_energy'] = free_energy

    #assume that the free energy is always the best to use
    mol.energy = free_energy


def geometry(mol,line_number,line,output_lines,calculator):
    global step
    step += 1
    firstxyzline = line_number + 5
    coords = []
    for parseline_number in range(0,mol.natoms):
        center,atomic_number,atomic_type,x,y,z = output_lines[parseline_number + firstxyzline].split()
        coords.append([float(x),float(y),float(z)])

    coords = np.array(coords)
    #update the mol object coordinates
    mol.coords = coords
    #propogate the coords to the xyz, xyzstring 
    mol.update_geometry()

    #add the geometry to the list of geometries
    if 'optimization_xyzs' in mol.properties:
        mol.properties['optimization_xyzs'].append(mol.ToXYZ(title='step{0}'.format(step),coords=coords,write=False))
    else:
        mol.properties['optimization_xyzs'] = [mol.ToXYZ(title='step{0}'.format(step),coords=coords,write=False)]

def non_opt(mol,line_number,line,output_lines,calculator):
    warning('''Non-Optimized Parameters - Gaussian did not converge the geometry for {0}

    The mol.warnings = ['Non-Optimized']
    
    If this is an opt_freq job, it is reccomended to submit again using
    the 'oldchk' keyword with this calculation's .chk file and
    reading the force constants opt='readfc'
    
    This calculation's .chk is {0}/{0}.chk
    
    '''.format(calculator.jobname))

    mol.warnings.append('Non-Optimized')

def scf_error(mol,line_number,line,output_lines,calculator):
    warning('''SCF ERROR - Gaussian did not converge the energy for {0}

    The mol.warnings = ['SCF_Error']
    
    submitting with scf=qc is reccomended 
    if you are confident in the input structure
    
    '''.format(calculator.jobname))

    mol.warnings.append('SCF_Error')

def frequencies(mol,line_number,line,output_lines,calculator):
    first_freq_line = line_number + 6
    freq_line = output_lines[first_freq_line].split()[2:]
    freqs = [float(freq) for freq in freq_line]

    #need to step forward in the file grabbing the frequencies
    #from every 7 + natoms lines
    freq_line_number = first_freq_line + 7 + mol.natoms
    freq_line = output_lines[first_freq_line]
    while re.search('Frequencies', freq_line):
        freq_line = output_lines[first_freq_line].split()[2:]
        freqs.extend([float(freq) for freq in freq_line])
        freq_line_number += 7
        freq_line_number += mol.natoms
        freq_line = output_lines[freq_line_number]
    
    freqs = np.array(freqs)

    mol.properties['frequencies'] = freqs

    if freqs[1] < 0:
        warning('''Saddle point found for {0}

        The mol.warnings = ['negative_frequency']
        
        It's reccomended to resubmit reading the force constants in
        
        '''.format(calculator.jobname))

        mol.warnings.append('negative_frequency')

    elif freqs[0] < 0:
        warning('''Negative frequency present for {0}

        The mol.warnings = ['negative_frequency']
        
        Ignore if this is a TS, otherwise it's reccomended
        to resubmit reading the force constants in
        
        '''.format(calculator.jobname))

        mol.warnings.append('negative_frequency')

def tddft_energy(mol,line_number,line,output_lines,calculator):
    mol.energy = float(line.split()[-1])
    mol.properties['tddft_energy'] = float(line.split()[-1])


def opt(outputfile):
    '''Properties to find
            1. electornic energy
            2. thermal energies
            3. coordinates
            4. homo
            5. lumo
            6. homo -1
            7. lumo +1
            8. frequencies
            9. list of optimization energies
            10. frequency vectors
            11. optimization geometries
                        1. electronic energy
            2. coordinates
            3. homo
            4. lumo
            5. homo -1
            6. lumo +1
            7. list of irc energies
            8. irc geometries
            Properties to find
        1. electornic energy
        2. homo
        3. lumo
        4. homo-1
        5. lumo +1'''
    pass
            



