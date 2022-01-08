#!/usr/bin/env python
#SBATCH --job-name=testmol_conf_search
#SBATCH --output=out
#SBATCH --error=error
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=23:59:59

###################################################################
#Start of the workflow
from workflowV2 import molecule
from workflowV2.calculator import Run
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2.software.CREST import CREST
from workflowV2 import message
import pandas as pd
import numpy as np
from datetime import datetime
from argparse import ArgumentParser
import os

###################################################################
def setup(mol,name, constraints=[],charge=0,mult=1,logfile='workflow.log'):
    """Get things prepared before starting the workflow"""

    #write output to log file
    message.logtofile(logfile)

    #look for checkpointed workflow
    if os.path.exists(name+'.chk'):
        mol = molecule.ReadCheckpoint(mol.replace('.xyz','.chk'))
        message.log('''##################################################################\nFound checkpoined mol for restart at step {1} - {0}\n##################################################################
            '''.format(datetime.now(), mol.tags['step']), time=False)


    else:
        mol = molecule.XYZToMol(mol,tags={'step':0},charge=charge,mult=mult,constraints=constraints)
        message.log('''##################################################################\nStarting workflow at - {0}\n##################################################################
        '''.format(datetime.now()), time=False)

        constraint_type = {2:'bond',3:'angle',4:'dihedral'}
        constraint_string  = ''.join(['\n    ' + constraint_type[len(constraint)] + ': ' + ','.join([str(c) for c in constraint]) for constraint in constraints ])
        message.log(f'Applying the following atom constraints:\n{constraint_string}\n\n',time=False)


    return(mol)

###################################################################

def conf_search(mol,method='B3LYP/6-31G',solvent=None,CREST_solvent=False,nproc=16,mem=120,partition='short',time='1-00:00',constraints=[],name='mol',max_confs=500,charge=0,mult=1,logfile='workflow.log',TS=False,n_final_opts=10):
    """Run the conformational search and optimization on an input TS xyz file"""

    ################################################################
    # Run the setup to create/read mol
    mol = setup(mol,name,constraints,charge,mult,logfile)
    chk = name+'.chk'
    csv = name+'_energies.csv'


    ################################################################
    # Run conformational search
    if mol.tags['step'] < 1:
        os.makedirs('CREST_conf_search',exist_ok=True)
        os.chdir('CREST_conf_search')

        use_solvent = solvent  #default to the regular solvent key
        if CREST_solvent:      #unless the CREST_solvent is given
            if CREST_solvent.lower() == 'none':   #if the CREST_solvent is 'none', then remove solvent
                use_solvent = False
            else:
                use_solvent = CREST_solvent      #otherwise set solvent to CREST_solvent

        if use_solvent:
            conformer_calculator = CREST(mol,jobname=f'{name}-confsearch',runtype='confsearch',
                                     gbsa=use_solvent,max_confs=max_confs,
                                     nproc=nproc,mem=mem,partition=partition,time=time)
        else:
            conformer_calculator = CREST(mol,jobname=f'{name}-confsearch',runtype='confsearch',
                                     max_confs=max_confs,
                                     nproc=nproc,mem=mem,partition=partition,time=time)


        mol = Run(conformer_calculator)

        message.log(f'{len(mol.conformers)} conformers found by CREST')

        #checkpoint the workflow
        os.chdir('../')

        #Write out a table of conformer energies
        conf_energy_pairs = [[idx,conf.energy] for idx,conf in enumerate(mol.conformers)]
        df = pd.DataFrame(conf_energy_pairs,columns=['CREST_Conformer','CREST_energy'])
        df.to_csv(csv,index=False)

        mol.tags['step'] = 1
        mol.checkpoint(chk)


    ################################################################
    # Remove the constraints and re-rank conformers
    if mol.tags['step'] < 2:
        os.makedirs('SP_conformer_ranking',exist_ok=True)
        os.chdir('SP_conformer_ranking')

        #Remove constraints
        mol.constraints = []
        for idx,conf in enumerate(mol.conformers):
            conf.constraints = []
            conf.tags['CREST_Conformer'] = idx   #add an idx tag to track where it ends up after reshuffling

        if solvent is not None:
            mol.RefineConformers(GAUSSIAN, jobname=f'{name}-sp',runtype='sp',
                             method=method,scrf=f'iefpcm,solvent={solvent}',
                             nproc=nproc,mem=mem,partition=partition,time=time)

        else:
            mol.RefineConformers(GAUSSIAN, jobname=f'{name}-sp',runtype='sp',
                             method=method,
                             nproc=nproc,mem=mem,partition=partition,time=time)


        #Get mapping of CREST_conformer -> re-ranked_conformer
        conf_mapping = {conf.tags['CREST_Conformer']:idx for idx,conf in enumerate(mol.conformers)}

        for idx,conf in enumerate(mol.conformers):
            conf.tags['Re-Ranked_Conformer'] = idx   #add an idx tag to track where it ends up after reshuffling

        #checkpoint the workflow
        os.chdir('../')

        #Write re-ranked energies to table
        df = pd.read_csv(csv)
        df['Re-Ranked_Conformer'] = df['CREST_Conformer'].map(conf_mapping)
        df.sort_values(['Re-Ranked_Conformer','CREST_Conformer'],inplace=True)
        df['Re-Ranked_SP_Energy'] = [conf.properties['electronic_energy'] for conf in mol.conformers]
        df.to_csv(csv,index=False)

        mol.tags['step'] = 2
        mol.checkpoint(chk)


    ################################################################
    # Optimize the lowest 10 conformers

    if mol.tags['step'] < 3:
        os.makedirs('final_optimizations', exist_ok=True)
        os.chdir('final_optimizations')


        #Cut conformer list to just the lowest 10
        if len(mol.conformers) > n_final_opts:
            mol.conformers = mol.conformers[0:n_final_opts]


        if TS:
            opt_arguments = 'ts, noeigen, recalcfc=10, maxstep=10'
        else:
            opt_arguments = 'recalcfc=15'

        if solvent is not None:
            mol.RefineConformers(GAUSSIAN,jobname=f'{name}-opt',runtype='opt_freq',
                             method=method,opt=opt_arguments,TS=TS,scrf=f'iefpcm,solvent={solvent}',
                             tries=3,nproc=nproc,mem=mem,partition=partition,time=time)
        else:
            mol.RefineConformers(GAUSSIAN,jobname=f'{name}-opt',runtype='opt_freq',
                             method=method,opt=opt_arguments,TS=TS,
                             tries=3,nproc=nproc,mem=mem,partition=partition,time=time)



        #Get mapping of CREST_conformer -> re-ranked_conformer
        conf_mapping = {conf.tags['Re-Ranked_Conformer']:idx for idx,conf in enumerate(mol.conformers)}

        # checkpoint the workflow
        os.chdir('../')

        #Write the final energies to table
        df = pd.read_csv(csv)
        df['Final_Conformer'] = df['Re-Ranked_Conformer'].map(conf_mapping)
        df.sort_values(['Final_Conformer','Re-Ranked_Conformer','CREST_Conformer'],inplace=True,na_position='last') #put null values last
        df['Final_Free_Energy'] = [mol.conformers[idx].properties['free_energy'] if idx < n_final_opts else np.nan for idx in range(len(df))]
        df['Final_Relative_Energy_kcal'] = [mol.conformers[idx].energy if idx < n_final_opts else np.nan for idx in range(len(df))]
        df.to_csv(csv,index=False)

        mol.tags['step'] = 3
        mol.checkpoint(name + '.chk')


    ################################################################
    # Write output
    if mol.tags['step'] < 4:
        mol.conformers[0].ToXYZ(f'{name}-lowestconf.xyz')
        mol.ConformersToXYZ(f'{name}-lowest{n_final_opts}confs.xyz')


##################################################################

##################################################################

##################################################################

if __name__ == '__main__':

    #parse arguments
    parser = ArgumentParser('Generalized Conformational Search Routine',fromfile_prefix_chars='@')
    parser.add_argument('--xyz',type=str,dest='mol',required=True,help='molecule input file in xyz format')
    parser.add_argument('-o','--logfile',type=str,dest='logfile',default='guess',help='name of file to log workflow output, defaults to guessing a name based on the xyz file')
    parser.add_argument('-n','--name',type=str,dest='name',default='guess',help='molecule name for calculation input/output files etc, defaults to guessing a name based on the xyz file')
    parser.add_argument('--maxconfs',dest='max_confs',type=int,default=500,help='maximim number of conformers to retrieve from CREST')
    parser.add_argument('-q','--charge',dest='charge',type=int,default=0,help='molecule charge')
    parser.add_argument('-m','--mult',dest='mult',type=int,default=1,help='molecule multiplicity')
    parser.add_argument('-b','--bond',dest='bonds',nargs=2,action='append',default=[],help='2 atom indicies to constrain a distance, 0 indexed - only applied to CREST conformational search, not subsequent optimizations')
    parser.add_argument('-a','--angle', dest='angles',nargs=3,action='append',default=[],help='3 atom indicies to constain an angle, 0 indexed - only applied to CREST conformational search, not subsequent optimizations')
    parser.add_argument('-d','--dihedral',dest='dihedrals',nargs=4,action='append',default=[],help='4 atom indicies to constrain a dihedral, 0 indexed - only applied to CREST conformational search, not subsequent optimizations')
    parser.add_argument('--TS',dest='TS',action='store_true',default=False,help='Flag to treat the input as a transition state structure')
    parser.add_argument('--final_opts',dest='n_final_opts',type=int,default=10,help='Maximum number of conformers to keep for the final full optimization')
    parser.add_argument('--solvent',dest='solvent', type=str,default=None,help='Solvent to use for calculations, if not specified will run in gas phase. See https://gaussian.com/scrf/ and https://xtb-docs.readthedocs.io/en/latest/gbsa.html?highlight=solvent for availible solvents')
    parser.add_argument('--CREST_solvent',dest='CREST_solvent', type=str,default=False,help='Separate solvent parameter to use for CREST conformational search only - can be set to "None". To be used if the final desired solvent is not availible in CREST, for example. See https://xtb-docs.readthedocs.io/en/latest/gbsa.html?highlight=solvent for availible solvents.')
    parser.add_argument('--method',dest='method',type=str,default='b3lyp/6-31G',help='Calculation method (dft functional and basis set for example)')
    parser.add_argument('--nproc',dest='nproc',type=int,default=16,help='Number of processors to use for calculations')
    parser.add_argument('--mem',dest='mem',type=int,default=120,help='Memory to be allocated to each calculation (Gb)')
    parser.add_argument('--partition',dest='partition',type=str,default='short',help='Partition(s) to run calculations on')
    parser.add_argument('--time',dest='time',type=str,default='1-00:00',help='Time limit for each individual calculation')


    arguments = parser.parse_args()

    #name the log file after the xyz file
    if arguments.logfile == 'guess':
        arguments.logfile = arguments.mol.replace('.xyz','.log')
        if arguments.logfile == arguments.mol:
            arguments.logfile = 'workflow.log'

    #name the name after the xyz file
    if arguments.name == 'guess':
        arguments.name = arguments.mol.split('.')[0]

    #ignore any empty constraint lists
    constraints = []
    for constraint in [arguments.bonds,arguments.angles,arguments.dihedrals]:
        if len(constraint) > 0:
            for const in constraint:
                constraints.append([int(c.strip()) for c in const])

    conf_search(mol=arguments.mol,
                solvent=arguments.solvent,
                CREST_solvent=arguments.CREST_solvent,
                method=arguments.method,
                nproc=arguments.nproc,
                mem=arguments.mem,
                partition=arguments.partition,
                time=arguments.time,
                constraints=constraints,
                name=arguments.name,
                max_confs=arguments.max_confs,
                charge=arguments.charge,
                mult=arguments.mult,
                logfile=arguments.logfile,
                TS=arguments.TS,
                n_final_opts=arguments.n_final_opts)

    message.log('Done!')

