from workflowV2 import molecule
from workflowV2.calculator import RunBatch
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2 import message
import pandas as pd
import os
from datetime import datetime
from .util import open_csv,unlock_csv,look_for_restart

#function to run S1 reactant optimizations

def S1_reactant_flow(mols,mol_names,output_energies,solvent,dftmethod,preopt_method,tddftmethod,nproc,mem,time,partition,mol_name):

    #now that we are in an independent branch, we need to locally checkpoint progress
    #use the 'S1s' mol objects to checkpoint this branch

    if os.path.exists(mol_names[0] + '-S1-reactant.chk'):
        S1s = [molecule.ReadCheckpoint(mol + '-S1-reactant.chk') for mol in mol_names]
        message.log('''##################################################################\nfound checkpoined S1s for restart - {0}\n##################################################################

    '''.format(datetime.now()),time=False)

    else:
        #if no checkpoint present, create the S1s from the mols
        #change the charge and multiplicity for each
        S1s = [mol.copy() for mol in mols]

        for S1 in S1s:
            S1.tags['step'] = 0

##################################################################

    if S1s[0].tags['step'] < 1:

        os.makedirs('S1-reactant-preopt',exist_ok=True)
        os.chdir('S1-reactant-preopt')

        calculators = [GAUSSIAN(S1,runtype='opt',jobname=name+'-S1-reactant-preopt',
                                    method=preopt_method,scrf=f'iefpcm,solvent={solvent}',
                                    opt='MaxCycles=15',
                                    td='Root=1,NStates=3',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    partition=partition) for index,S1,name in zip(range(len(S1s)),S1s,mol_names)]

        S1s = RunBatch(calculators,tries=1,jobname=f'{mol_name}-S1-reactant-preopts',ignore=True)

        for S1,name in zip(S1s,mol_names):
            S1.checkpoint('{0}-S1_reactant-preopt.mol'.format(name))

        os.chdir('../')

        for S1,name in zip(S1s,mol_names):
            S1.tags['step'] = 1
            S1.checkpoint('{0}-S1-reactant.chk'.format(name))

##################################################################

    if S1s[0].tags['step'] < 2:

        os.makedirs('S1-reactant',exist_ok=True)
        os.chdir('S1-reactant')

        oldchk,geom,guess = look_for_restart(mol_names,'S1-reactant')
        calculators = [GAUSSIAN(S1,runtype='opt_freq',jobname=name+'-S1-reactant',
                                    method=dftmethod,scrf=f'iefpcm,solvent={solvent}',
                                    td='Root=1,NStates=3',
                                    opt='recalcfc=10',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    oldchk=oldchk[index],
                                    geom=geom[index],
                                    guess=guess[index],
                                    partition=partition) for index,S1,name in zip(range(len(S1s)),S1s,mol_names)]

        S1s = RunBatch(calculators,tries=3,jobname=f'{mol_name}-S1-reactants')

        for S1,name in zip(S1s,mol_names):
            S1.checkpoint('{0}-S1_reactant.mol'.format(name))

        os.chdir('../')

        df = open_csv(output_energies,'S1-reactant')
        df['S1_reactant_free_energy'] = [S1.energy for S1 in S1s]

        #if the anions and cations are done, calculate the redox potentials
        if 'S1_reduced_free_energy' in df.columns:
             df['S1_reduction_free_energy (eV)'] = (df['S1_reduced_free_energy'] - df['S1_reactant_free_energy']) * 27.211
             df['S1_reduction_potential (V)'] = - df['S1_reduction_free_energy'] - 4.44
        if 'S1_oxidized_free_energy' in df.columns:
             df['S1_oxidation_free_energy (eV)'] = (df['S1_oxidized_free_energy'] - df['S1_reactant_free_energy']) * 27.211
             df['S1_oxidation_potential (V)'] = - df['S1_oxidation_free_energy'] - 4.44

        df.to_csv(output_energies,index=False)
        unlock_csv(output_energies,'S1-reactant')

        for S1,name in zip(S1s,mol_names):
            S1.tags['step'] = 2
            S1.checkpoint('{0}-S1-reactant.chk'.format(name))

##################################################################

    message.log('S1 reactant sub-process complete')
    return('complete')


