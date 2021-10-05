from workflowV2 import molecule
from workflowV2.calculator import RunBatch
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2 import message
import pandas as pd
import os
from datetime import datetime
from .util import open_csv,unlock_csv,look_for_restart,par_proc,radical_electrons

#################################################################

def S0_reduction_flow(reduceds,mol_names,dftmethod,solvent,partition,output_energies,nproc,mem,time,mol_name):

    if os.path.exists(mol_names[0] + '-S0-reduced.chk'):
        reduceds = [molecule.ReadCheckpoint(reduced + '-S0-reduced.chk') for reduced in mol_names]
        message.log('''##################################################################\nfound checkpoined S0-reduceds for restart - {0}\n##################################################################
        '''.format(datetime.now()),time=False)

    else:
        reduceds = [reduced.copy() for reduced in reduceds]
        
        for reduced in reduceds:
            reduced.tags['step'] = 0

#################################################################

    if reduceds[0].tags['step'] < 1:

        os.makedirs('S0-reduced',exist_ok=True)
        os.chdir('S0-reduced')

        oldchk,geom,guess = look_for_restart(mol_names,'S0-reduced')
        calculators = [GAUSSIAN(reduced,runtype='opt_freq',jobname=name+'-S0-reduced',
                                    method=dftmethod,scrf=f'iefpcm,solvent={solvent}',
                                    opt='recalcfc=10',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    oldchk=oldchk[index],
                                    geom=geom[index],
                                    guess=guess[index],
                                    partition=partition) for index,reduced,name in zip(range(len(reduceds)),reduceds,mol_names)]

        reduceds = RunBatch(calculators,tries=3,jobname=f'{mol_name}-S0-reduceds')

        for reduced,name in zip(reduceds,mol_names):
            reduced.checkpoint('{0}-reduced.mol'.format(name))

        os.chdir('../')

        df = open_csv(output_energies,'S0-reduced')
        df['S0-reduced_free_energy'] = [reduced.energy for reduced in reduceds]
        df['S0_reduction_free_energy (eV)'] = (df['S0-reduced_free_energy'] - df['S0_reactant_free_energy']) * 27.211
        df['S0_reduction_potential (V)'] =  - df['S0_reduction_free_energy (eV)'] - 4.44
        df.to_csv(output_energies,index=False)
        unlock_csv(output_energies,'S0-reduced')

        for reduced,name in zip(reduceds,mol_names):
            reduced.tags['step'] = 1
            reduced.checkpoint('{0}-S0-reduced.chk'.format(name))

    message.log('S0 reduced sub-process complete')
    return('complete')


#################################################################

def S1_reduction_flow(reduceds,mol_names,dftmethod,preopt_method,solvent,partition,output_energies,nproc,mem,time,mol_name):

    if os.path.exists(mol_names[0] + '-S1-reduced.chk'):
        reduceds = [molecule.ReadCheckpoint(reduced + '-S1-reduced.chk') for reduced in mol_names]
        message.log('''##################################################################\nfound checkpoined S1-reduceds for restart - {0}\n##################################################################
        '''.format(datetime.now()),time=False)

    else:
        reduceds = [reduced.copy() for reduced in reduceds]
        
        for reduced in reduceds:
            reduced.tags['step'] = 0

#################################################################

    if reduceds[0].tags['step'] < 1:

        os.makedirs('S1-reduced-preopt',exist_ok=True)
        os.chdir('S1-reduced-preopt')

        oldchk,geom,guess = look_for_restart(mol_names,'S1-reduced-preopt')
        calculators = [GAUSSIAN(reduced,runtype='opt',jobname=name+'-S1-reduced-preopt',
                                    method=preopt_method,scrf=f'iefpcm,solvent={solvent}',
                                    opt='MaxCycles=15',
                                    td='Root=1,NStates=3',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    partition=partition) for index,reduced,name in zip(range(len(reduceds)),reduceds,mol_names)]

        reduceds = RunBatch(calculators,tries=1,jobname=f'{mol_name}-S1-reduced-preopts',ignore=True)

        for reduced,name in zip(reduceds,mol_names):
            reduced.checkpoint('{0}-S1-reduced-preopts.mol'.format(name))

        os.chdir('../')

        for reduced,name in zip(reduceds,mol_names):
            reduced.tags['step'] = 1
            reduced.checkpoint('{0}-S1-reduced.chk'.format(name))

#################################################################

    if reduceds[0].tags['step'] < 2:

        os.makedirs('S1-reduced',exist_ok=True)
        os.chdir('S1-reduced')

        oldchk,geom,guess = look_for_restart(mol_names,'S1-reduced')
        calculators = [GAUSSIAN(reduced,runtype='opt_freq',jobname=name+'-S1-reduced',
                                    method=dftmethod,scrf=f'iefpcm,solvent={solvent}',
                                    td='Root=1,NStates=3',
                                    opt='recalcfc=10',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    oldchk=oldchk[index],
                                    geom=geom[index],
                                    guess=guess[index],
                                    partition=partition) for index,reduced,name in zip(range(len(reduceds)),reduceds,mol_names)]

        reduceds = RunBatch(calculators,tries=3,jobname=f'{mol_name}-S1-reduceds')

        for reduced,name in zip(reduceds,mol_names):
            reduced.checkpoint('{0}-S1-reduced.mol'.format(name))

        os.chdir('../')

        df = open_csv(output_energies,'S1-reduced')
        df['S1_reduced_free_energy'] = [reduced.energy for reduced in reduceds]

        #if the neutral S1 is done, calculate the S1 reduced potential
        if 'S1_reactant_free_energy' in df.columns:
            df['S1_reduction_free_energy (eV)'] = (df['S1_reduced_free_energy'] - df['S1_reactant_free_energy']) * 27.211
            df['S1_reduction_potential (V)'] = - df['S1_reduction_free_energy (eV)'] - 4.44

        df.to_csv(output_energies,index=False)
        unlock_csv(output_energies,'S1-reduced')

        for reduced,name in zip(reduceds,mol_names):
            reduced.tags['step'] = 2
            reduced.checkpoint('{0}-S1-reduced.chk'.format(name))

    message.log('S1 reduced sub-process complete')
    return('complete')


#################################################################

#function to run reduced ground state and reduced S1 optimizations
def reduction_flow(mols,mol_names,output_energies,solvent,dftmethod,preopt_method,tddftmethod,nproc,mem,time,partition,mol_name):

    #now that we are in an independent branch, we need to locally checkpoint progress
    #use the 'reduced' mol objects to checkpoint this branch

    if os.path.exists(mol_names[0] + '-reduced.chk'):
        reduceds = [molecule.ReadCheckpoint(mol + '-reduced.chk') for mol in mol_names]
        message.log('''##################################################################\nfound checkpoined reduceds for restart - {0}\n##################################################################

    '''.format(datetime.now()),time=False)

    else:
        #if no checkpoint present, create the reduceds from the mols
        #change the charge and multiplicity for each
        reduceds = [mol.copy() for mol in mols]

        for reduced in reduceds:
            reduced.charge -= 1

            n_radical_electrons = radical_electrons(reduced.rdkitmol)
            n_radical_electrons += 1
            multiplicity = abs(n_radical_electrons) + 1


            reduced.mult = multiplicity #int(((oxidized.charge / 2) * 2) + 1)

            reduced.tags['step'] = 0

##################################################################

    if reduceds[0].tags['step'] < 1:

        os.makedirs('reduced-preopt',exist_ok=True)
        os.chdir('reduced-preopt')

        calculators = [GAUSSIAN(reduced,runtype='opt',jobname=name+'-reduced-preopt',
                                    method=preopt_method,scrf=f'iefpcm,solvent={solvent}',
                                    opt='MaxCycles=15',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    partition=partition) for index,reduced,name in zip(range(len(reduceds)),reduceds,mol_names)]

        reduceds = RunBatch(calculators,tries=1,jobname=f'{mol_name}-reduced-preopts',ignore=True)

        for reduced,name in zip(reduceds,mol_names):
            reduced.checkpoint('{0}-reduced-preopt.mol'.format(name))

        os.chdir('../')

        for reduced,name in zip(reduceds,mol_names):
            reduced.tags['step'] = 1
            reduced.checkpoint('{0}-reduced.chk'.format(name))

##################################################################

# NOW THAT WE HAVE A CLEANED UP ANION GEOMETRY,
# LETS SPLIT OFF PROCESSES FOR GROUND AND EXCITED STATES

    if reduceds[0].tags['step'] < 2:

        task_list = [

            {'func':S0_reduction_flow,'tasks':[{'reduceds':reduceds,'mol_names':mol_names,'output_energies':output_energies,
                                                        'solvent':solvent,'dftmethod':dftmethod,
                                                        'nproc':nproc,'mem':mem,'time':time,'partition':partition,'mol_name':mol_name}]},

            {'func':S1_reduction_flow,'tasks':[{'reduceds':reduceds,'mol_names':mol_names,'output_energies':output_energies,
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
            message.warning(f'Failed calculations in the {failed_branches} reduction-sub-processes!!! ... exiting')
            raise ValueError('sub-process did not return complete')
  
        for reduced,name in zip(reduceds,mol_names):
            reduced.tags['step'] = 2
            reduced.checkpoint('{0}-reduced.chk'.format(name)) 

##################################################################

    message.log('Anion sub-process complete')
    return('complete')


