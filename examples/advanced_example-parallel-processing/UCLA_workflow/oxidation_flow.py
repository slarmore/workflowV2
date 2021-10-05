from workflowV2 import molecule
from workflowV2.calculator import RunBatch
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2 import message
import pandas as pd
import os
from datetime import datetime
from .util import open_csv,unlock_csv,look_for_restart,par_proc,radical_electrons

#################################################################

def S0_oxidation_flow(oxidizeds,mol_names,dftmethod,solvent,partition,output_energies,nproc,mem,time,mol_name):

    if os.path.exists(mol_names[0] + '-S0-oxidized.chk'):
        oxidizeds = [molecule.ReadCheckpoint(oxidized + '-S0-oxidized.chk') for oxidized in mol_names]
        message.log('''##################################################################\nfound checkpoined S0-oxidizeds for restart - {0}\n##################################################################
        '''.format(datetime.now()),time=False)

    else:
        oxidizeds = [oxidized.copy() for oxidized in oxidizeds]
        
        for oxidized in oxidizeds:
            oxidized.tags['step'] = 0

#################################################################

    if oxidizeds[0].tags['step'] < 1:

        os.makedirs('S0-oxidized',exist_ok=True)
        os.chdir('S0-oxidized')

        oldchk,geom,guess = look_for_restart(mol_names,'S0-oxidized')
        calculators = [GAUSSIAN(oxidized,runtype='opt_freq',jobname=name+'-S0-oxidized',
                                    method=dftmethod,scrf=f'iefpcm,solvent={solvent}',
                                    opt='recalcfc=10',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    oldchk=oldchk[index],
                                    geom=geom[index],
                                    guess=guess[index],
                                    partition=partition) for index,oxidized,name in zip(range(len(oxidizeds)),oxidizeds,mol_names)]

        oxidizeds = RunBatch(calculators,tries=3,jobname=f'{mol_name}-S0-oxidizeds')

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.checkpoint('{0}-oxidized.mol'.format(name))

        os.chdir('../')

        df = open_csv(output_energies,'S0-oxidized')
        df['S0_oxidized_free_energy'] = [oxidized.energy for oxidized in oxidizeds]
        df['S0_oxidation_free_energy (eV)'] = (df['S0_reactant_free_energy'] - df['S0_oxidized_free_energy']) * 27.211
        df['S0_oxidation_potential (V)'] = - df['S0_oxidation_free_energy (eV)'] - 4.44
        df.to_csv(output_energies,index=False)
        unlock_csv(output_energies,'S0-oxidized')

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.tags['step'] = 1
            oxidized.checkpoint('{0}-S0-oxidized.chk'.format(name))

    message.log('S0 oxidized sub-process complete')
    return('complete')


#################################################################

def S1_oxidation_flow(oxidizeds,mol_names,dftmethod,preopt_method,solvent,partition,output_energies,nproc,mem,time,mol_name):

    if os.path.exists(mol_names[0] + '-S1-oxidized.chk'):
        oxidizeds = [molecule.ReadCheckpoint(oxidized + '-S1-oxidized.chk') for oxidized in mol_names]
        message.log('''##################################################################\nfound checkpoined S1-oxidizeds for restart - {0}\n##################################################################
        '''.format(datetime.now()),time=False)

    else:
        oxidizeds = [oxidized.copy() for oxidized in oxidizeds]
        
        for oxidized in oxidizeds:
            oxidized.tags['step'] = 0

#################################################################

    if oxidizeds[0].tags['step'] < 1:

        os.makedirs('S1-oxidized-preopt',exist_ok=True)
        os.chdir('S1-oxidized-preopt')

        oldchk,geom,guess = look_for_restart(mol_names,'S1-oxidized-preopt')
        calculators = [GAUSSIAN(oxidized,runtype='opt',jobname=name+'-S1-oxidized-preopt',
                                    method=preopt_method,scrf=f'iefpcm,solvent={solvent}',
                                    opt='MaxCycles=15',
                                    td='Root=1,NStates=3',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    partition=partition) for index,oxidized,name in zip(range(len(oxidizeds)),oxidizeds,mol_names)]

        oxidizeds = RunBatch(calculators,tries=1,jobname=f'{mol_name}-S1-oxidized-preopts',ignore=True)

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.checkpoint('{0}-S1_oxidized-preopt.mol'.format(name))

        os.chdir('../')

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.tags['step'] = 1
            oxidized.checkpoint('{0}-S1-oxidized.chk'.format(name))

#################################################################

    if oxidizeds[0].tags['step'] < 2:

        os.makedirs('S1-oxidized',exist_ok=True)
        os.chdir('S1-oxidized')

        oldchk,geom,guess = look_for_restart(mol_names,'S1-oxidized')
        calculators = [GAUSSIAN(oxidized,runtype='opt_freq',jobname=name+'-S1-oxidized',
                                    method=dftmethod,scrf=f'iefpcm,solvent={solvent}',
                                    td='Root=1,NStates=3',
                                    opt='recalcfc=10',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    oldchk=oldchk[index],
                                    geom=geom[index],
                                    guess=guess[index],
                                    partition=partition) for index,oxidized,name in zip(range(len(oxidizeds)),oxidizeds,mol_names)]

        oxidizeds = RunBatch(calculators,tries=3,jobname=f'{mol_name}-S1-oxidizeds')

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.checkpoint('{0}-S1_oxidized.mol'.format(name))

        os.chdir('../')

        df = open_csv(output_energies,'S1-oxidized')
        df['S1_oxidized_free_energy'] = [oxidized.energy for oxidized in oxidizeds]

        #if the neutral S1 is done, calculate the S1 oxidized potential
        if 'S1_reactant_free_energy' in df.columns:
            df['S1_oxidation_free_energy (eV)'] = (df['S1_reactant_free_energy'] - df['S1_oxidized_free_energy']) * 27.211
            df['S1_oxidation_potential (V)'] = - df['S1_oxidation_free_energy (eV)'] - 4.44

        df.to_csv(output_energies,index=False)
        unlock_csv(output_energies,'S1-oxidized')

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.tags['step'] = 2
            oxidized.checkpoint('{0}-S1-oxidized.chk'.format(name))

    message.log('S1 oxidized sub-process complete')
    return('complete')


#################################################################

#function to run oxidized ground state and oxidized S1 optimizations
def oxidation_flow(mols,mol_names,output_energies,solvent,dftmethod,preopt_method,tddftmethod,nproc,mem,time,partition,mol_name):

    #now that we are in an independent branch, we need to locally checkpoint progress
    #use the 'oxidized' mol objects to checkpoint this branch

    if os.path.exists(mol_names[0] + '-oxidized.chk'):
        oxidizeds = [molecule.ReadCheckpoint(mol + '-oxidized.chk') for mol in mol_names]
        message.log('''##################################################################\nfound checkpoined oxidizeds for restart - {0}\n##################################################################

    '''.format(datetime.now()),time=False)

    else:
        #if no checkpoint present, create the oxidizeds from the mols
        #change the charge and multiplicity for each
        oxidizeds = [mol.copy() for mol in mols]

        for oxidized in oxidizeds:
            oxidized.charge += 1

            n_radical_electrons = radical_electrons(oxidized.rdkitmol)
            n_radical_electrons -= 1
            multiplicity = abs(n_radical_electrons) + 1


            oxidized.mult = multiplicity #int(((oxidized.charge / 2) * 2) + 1)
            oxidized.tags['step'] = 0

##################################################################

    if oxidizeds[0].tags['step'] < 1:

        os.makedirs('oxidized-preopt',exist_ok=True)
        os.chdir('oxidized-preopt')

        calculators = [GAUSSIAN(oxidized,runtype='opt',jobname=name+'-oxidized-preopt',
                                    method=preopt_method,scrf=f'iefpcm,solvent={solvent}',
                                    opt='MaxCycles=15',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    partition=partition) for index,oxidized,name in zip(range(len(oxidizeds)),oxidizeds,mol_names)]

        oxidizeds = RunBatch(calculators,tries=1,jobname=f'{mol_name}-oxidized-preopts',ignore=True)

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.checkpoint('{0}-oxidized-preopt.mol'.format(name))

        os.chdir('../')

        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.tags['step'] = 1
            oxidized.checkpoint('{0}-oxidized.chk'.format(name))

##################################################################

# NOW THAT WE HAVE A CLEANED UP ANION GEOMETRY,
# LETS SPLIT OFF PROCESSES FOR GROUND AND EXCITED STATES

    if oxidizeds[0].tags['step'] < 2:

        task_list = [

            {'func':S0_oxidation_flow,'tasks':[{'oxidizeds':oxidizeds,'mol_names':mol_names,'output_energies':output_energies,
                                                        'solvent':solvent,'dftmethod':dftmethod,
                                                        'nproc':nproc,'mem':mem,'time':time,'partition':partition,'mol_name':mol_name}]},

            {'func':S1_oxidation_flow,'tasks':[{'oxidizeds':oxidizeds,'mol_names':mol_names,'output_energies':output_energies,
                                                        'solvent':solvent,'dftmethod':dftmethod,'preopt_method':preopt_method,
                                                        'nproc':nproc,'mem':mem,'time':time,'partition':partition,'mol_name':mol_name}]}
                    ]


        results = par_proc(task_list,num_cpus=2)
        failed_branches = []
        for result in results:
            for func,status in result.items():
                if status != 'complete':
                    failed_branches.append(func)
        if len(failed_branches) > 0:
            message.warning(f'Failed calculations in the {failed_branches} oxidation-sub-processes!!! ... exiting')
            raise ValueError('sub-process did not return complete')
  
        for oxidized,name in zip(oxidizeds,mol_names):
            oxidized.tags['step'] = 2
            oxidized.checkpoint('{0}-oxidized.chk'.format(name)) 

##################################################################

    message.log('Oxidation sub-process complete')
    return('complete')


