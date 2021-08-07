#!/usr/bin/env python
#SBATCH --job-name=binol_workflow
#SBATCH --output=out
#SBATCH --error=error
#SBATCH --partition=lopez
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=5-00:00:00


solvent = 'toluene'
semiempericalmethod = 'pm7'
dftmethod = 'B97D/6-31G(d)'
nproc = 32
mem = 250


###################################################################

from workflowV2 import molecule
from workflowV2.calculator import Run
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2.software.CREST import CREST
from workflowV2 import message
import pickle
import pandas as pd
import sys
import os


input_file = sys.argv[1]

input_name = input_file.split('.')[0]

#log the output to a file
message.logtofile('{0}.log'.format(input_name))

##################################################################

#set up the molecule object
mol = molecule.XYZToMol(input_file,charge=-1)

#conformational search
os.makedirs('conf_search',exist_ok=True)
os.chdir('conf_search')

conformer_calculator = CREST(mol,jobname='{0}-confsearch'.format(input_name),runtype='confsearch',nproc=nproc,mem=mem,gbsa=solvent)

mol = Run(conformer_calculator)

os.chdir('../')

mol.tags['step'] = 1
with open('{0}.chk'.format(input_name),'wb') as chk:
    pickle.dump(mol,chk)
with open('{0}-confsearch.mol'.format(input_name),'wb') as chk:
    pickle.dump(mol,chk)


##################################################################
#Single points with TDDFT
os.makedirs('TDDFT-SP',exist_ok=True)
os.chdir('TDDFT-SP')

mol.RefineConformers(GAUSSIAN,jobname='{0}-TDDFT-sp'.format(input_name),
                             runtype='sp',method=method,td='Root=1,NStates=3',
                             scrf='iefpcm,solvent={0}'.format(solvent),
                             nproc=nproc,mem=mem)

os.chdir('../')

mol.tags['step'] = 2
with open('{0}.chk'.format(input_name),'wb') as chk:
    pickle.dump(mol,chk)
with open('{0}-TDDFT-SP.mol'.format(input_name),'wb') as chk:
    pickle.dump(mol,chk)



##################################################################
#optimize the lowest energy TDDFT conformer with ZIndo
os.makedirs('semiempirical',exist_ok=True)
os.chdir('semiempierical')

#just grab the lowest conformer
conf = mol.conformers[0]

conf.ToXYZ('{0}-starting_conf.xyz'.format(input_name))

calc = GAUSSIAN(conf,runtype='opt_freq',jobname='{0}-ZIndo'.format(input_name),
                method=semiempiricalmethod,
                scrf='iefpcm,solvent={0}'.format(solvent),
                nproc=nproc,mem=mem)

semiempirical = Run(calc)

semiempirical_energy = semiempirical.energy

ZIndo.ToXYZ('{0}-semiempirical.xyz'.format(input_name))

os.chdir('../')

mol.tags['step'] = 3
with open('{0}.chk','wb') as chk:
    pickle.dump(ZIndo,chk)

with open('{0}-semiempirical.mol','wb') as chk:
    pickle.dump(semiempirical,chk)


##################################################################
#optimize the lowest energy TDDFT conformer with ZIndo
os.makedirs('dft',exist_ok=True)
os.chdir('dft')

#just grab the lowest conformer
calc = GAUSSIAN(semiempirical,runtype='opt_freq',jobname='{0}-dft'.format(input_name),
                method=dftmethod,td='Root=1,NStates=3',
                scrf='iefpcm,solvent={0}'.format(solvent),
                nproc=nproc,mem=mem)

dft = Run(calc)

dft_energy = dft.energy

dft.ToXYZ('{0}-dft.xyz'.format(input_name))


os.chdir('../')

mol.tags['step'] = 3
with open('{0}.chk'.format(input_name),'wb') as chk:
    pickle.dump(dft,chk)

with open('{0}-dft.mol'.format(input_name),'wb') as chk:
    pickle.dump(dft,chk)


################################################################
#Write the output
energies = [[input_name,semiempirical_energy,dft_energy]]
energy_output = pd.DataFrame(energies,columns=['mol','semiempirical_energy','dft_energy'])
energy_output.to_csv('{0}-output_energies.csv'.format(input_name),index=False)

print('Done!')
