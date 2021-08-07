#!/usr/bin/env python
#SBATCH --job-name=conformational_search
#SBATCH --output=out
#SBATCH --error=error
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=5-00:00:00

########################################################################################################
#Settings for the user

#calculation specifications
solvent = 'water'
method = 'M062X/6-31+g(d,p)'

#workflow specifications
max_confs = 100                   #generate at most X conformers
max_optimize_confs = 10           #fully optimize the lowest X conformers for final ranking

#slurm options
nproc = 16
mem = 50
partition = 'short'
time = '1-00:00:00'

########################################################################################################
#Workflow description
 
              #######################################
              #       STEP 0: Input XYZ file        #
              #######################################

#                               |
#                               |  
#                               |
#                               V

              #######################################
              # STEP 1: CREST conformational search #
              #######################################
              # - Generate conformers with CREST

#                               |
#                               |  <max_conformers> on to next step
#                               |
#                               V

              #######################################
              # STEP 2: Rank conformers with DFT SP #
              #######################################
              # - For every conformer generated, 
              #   do a DFT SP with G16 to re-rank
              #   the conformers by energy

#                               |
#                               |  <max_optimize_conformers> on to next step
#                               |
#                               V

              #######################################
              # STEP 3: Rank conformers with DFT    #
              #         Optimized Free Energies     #
              #######################################
              # - For every conformer generated, 
              #   do a DFT SP with G16 to re-rank
              #   the conformers by energy

#                               |
#                               |  Lowest energy conformer
#                               |
#                               V

              #######################################
              #   Output: Lowest energy conformer   #
              #######################################
              # - Write XYZ of the lowest energy
              #   conformer
              # - Write xyz of all optimized conformers
              # - Write a csv of the conformer energies
    
    
########################################################################################################
#Start of script
        
from workflowV2 import molecule
from workflowV2.calculator import Run
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2.software.CREST import CREST
from workflowV2 import message
import pandas as pd
import sys
import os


#######################################
#       STEP 0: Input XYZ file        #
#######################################
input_file = sys.argv[1]
input_name = input_file.split('.')[0]

#log the output to a file
message.logtofile('{0}.log'.format(input_name))

#set up the molecule object
mol = molecule.XYZToMol(input_file)


#######################################
# STEP 1: CREST conformational search #
#######################################
os.makedirs('conf_search',exist_ok=True)
os.chdir('conf_search')

conformer_calculator = CREST(mol,max_confs=max_confs,
                             runtype='confsearch',gbsa=solvent,
                             nproc=nproc,mem=mem,partition=partition,time=time,
                             jobname='{0}-confsearch'.format(input_name))

mol = Run(conformer_calculator)

os.chdir('../')


#######################################
# STEP 2: Rank conformers with DFT SP #
#######################################
os.makedirs('SP',exist_ok=True)
os.chdir('SP')

mol.RefineConformers(GAUSSIAN,jobname='{0}-sp'.format(input_name),
                             runtype='sp',method=method,
                             scrf='iefpcm,solvent={0}'.format(solvent),
                             nproc=nproc,mem=mem,partition=partition,time=time)

os.chdir('../')


#######################################
# STEP 3: Rank conformers with DFT    #
#         Optimized Free Energies     #
#######################################
os.makedirs('OPT',exist_ok=True)
os.chdir('OPT')

#if there are more than max_optimize_conformers, then only keep max_optimize_conformers
if len(mol.conformers) > max_optimize_confs:
    mol.conformers = mol.conformers[0:max_optimize_confs]

mol.RefineConformers(GAUSSIAN,jobname='{0}-opt'.format(input_name),
                              runtype='opt_freq',method=method,
                              opt='recalcfc=10',
                              scrf='iefpcm,solvent={0}'.format(solvent),
                              nproc=nproc,mem=mem,partition=partition,time=time)

os.chdir('../')


#######################################
#   Output: Lowest energy conformer   #
#######################################
os.makedirs('lowest_energy_conformers',exist_ok=True)
os.chdir('lowest_energy_conformers')

#write the lowest energy xyz file
mol.conformers[0].ToXYZ('{0}-lowestconf.xyz'.format(input_name))

#write all of the optimized conformer xyzs
mol.ConformersToXYZ('{0}-allconf.xyz'.format(input_name))

#write out the absolute free energies and the relative energies
names_energies = [['{0}-conformer{1}'.format(input_name,index),conformer.properties['free_energy']] for index,conformer in enumerate(mol.conformers)]
energy_output = pd.DataFrame(names_energies,columns=['conformer','Free_energy'])
energy_output['relative_energy'] = [conformer.energy for conformer in mol.conformers]
energy_output.to_csv('{0}-output_energies.csv'.format(input_name),index=False)

os.chdir('../')

message.log('Workflow completed')
