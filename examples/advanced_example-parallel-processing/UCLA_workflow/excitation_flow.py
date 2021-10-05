from workflowV2 import molecule
from workflowV2.calculator import RunBatch
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2 import message
import pandas as pd
import os
from datetime import datetime
from .util import open_csv,unlock_csv,look_for_restart

#function to run tddft sp calculation for
#photophysics information


def excitation_flow(mols,mol_names,output_energies,solvent,dftmethod,tddftmethod,nproc,mem,time,partition,mol_name):

    #now that we are in an independent branch, we need to locally checkpoint progress
    #use the 'tddfts' mol objects to checkpoint this branch

    if os.path.exists(mol_names[0] + '-tddft.chk'):
        tddfts = [molecule.ReadCheckpoint(mol + '-tddft.chk') for mol in mol_names]
        message.log('''##################################################################\nfound checkpoined tddfts for restart - {0}\n##################################################################

    '''.format(datetime.now()),time=False)

    else:
        #if no checkpoint present, create the cations from the mols
        #change the charge and multiplicity for each
        tddfts = [mol.copy() for mol in mols]

        for tddft in tddfts:
            tddft.tags['step'] = 0

 ##################################################################

    if tddfts[0].tags['step'] < 1:

        os.makedirs('tddft',exist_ok=True)
        os.chdir('tddft')

        calculators = [GAUSSIAN(tddft,runtype='sp',jobname=name+'-tddft',
                                    method=tddftmethod,scrf=f'iefpcm,solvent={solvent}',
                                    td='NStates=10',
                                    nproc=nproc,
                                    mem=mem,
                                    time=time,
                                    partition=partition) for tddft,name in zip(tddfts,mol_names)]

        tddfts = RunBatch(calculators,tries=1,jobname='tddft')

        for tddft,name in zip(tddfts,mol_names):
            tddft.checkpoint('{0}-tddft.mol'.format(name))

        os.chdir('../')

        df = open_csv(output_energies,'tddft')

        #Need to unpack the dictionary of excited state properties into 
        #columns to add to our dataframe

        #crate 30 columns, 10 different states, energy, wavelength, oscillator strength per state
        for excited_state,state in enumerate(tddfts[0].properties['excited_states']):
            energies = []
            wavelengths = []
            oscillator_strengths = []
            for tddft in tddfts:
                energies.append(tddft.properties['excited_states'][excited_state]['energy'])
                wavelengths.append(tddft.properties['excited_states'][excited_state]['wavelength'])
                oscillator_strengths.append(tddft.properties['excited_states'][excited_state]['f'])
            df[f'S{excited_state + 1} energy (eV)'] = energies
            df[f'S{excited_state + 1}wavelength (nm)'] = wavelengths
            df[f'S{excited_state + 1} f'] = oscillator_strengths
        df.to_csv(output_energies,index=False)
        unlock_csv(output_energies,'tddft')
  
        for tddft,name in zip(tddfts,mol_names):
            tddft.tags['step'] = 1
            tddft.checkpoint('{0}-tddft.chk'.format(name))

##################################################################

    message.log('Excitation sub-process complete')
    return('complete')

