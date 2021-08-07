#!/usr/bin/env python
#SBATCH --job-name=hoop_workflow
#SBATCH --output=out-%A
#SBATCH --error=error-%A
#SBATCH --partition=lopez,short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=1-00:00:00

solvent_dict = {'crest':'DMSO','gaussian':'dimethylsulfoxide'}
dftmethod = 'M062X/6-31G'
nproc = 12
mem = 60

###################################################################

from workflowV2 import molecule
from workflowV2.calculator import Run as Run
from workflowV2.software.GAUSSIAN import GAUSSIAN as GAUSSIAN
from workflowV2.software.CREST import CREST
from workflowV2 import message
import pandas as pd
import sys
import os
import numpy as np


input_file = sys.argv[1]

input_name = input_file.split('.')[0]

#log the output to a file
message.logtofile('{0}.log'.format(input_name))

#write energies to file
output_energies = '{0}-output_energies.csv'.format(input_name)

#check if there was a checkpoint
mol = False
forward = False
reverse = False
product = False

if os.path.exists('{0}.chk'.format(input_name)):
    mol = molecule.ReadCheckpoint('{0}.chk'.format(input_name))
    message.log('found checkpoined mol')
    message.log('restarting at step {0}'.format(mol.tags['step']+1))

if not mol:
    atom1 = int(sys.argv[2])
    atom2 = int(sys.argv[3])
    atom3 = int(sys.argv[4])
    atom4 = int(sys.argv[5])

    mol = molecule.XYZToMol(input_file)
    mol.constraints = [[atom1,atom2],[atom3,atom4]]

    mol.tags['step'] = 0

##################################################################

if mol.tags['step'] < 1:
    #set up the molecule object
    #conformational search
    os.makedirs('conf_search',exist_ok=True)
    os.chdir('conf_search')

    conformer_calculator = CREST(mol,
                                 jobname='{0}-confsearch'.format(input_name),
                                 runtype='confsearch',
                                 nproc=nproc,
                                 mem=mem,
                                 gbsa=solvent_dict['crest'],
                                 max_confs=100)

    mol = Run(conformer_calculator)

    #release the constraints now that I have the conformers I needed
    mol.constraints = []
    for conformer in mol.conformers:
        conformer.constraints = []

    os.chdir('../')

    mol.tags['step'] = 1
    mol.checkpoint('{0}.chk'.format(input_name))
    mol.checkpoint('{0}-confsearch.mol'.format(input_name))

##################################################################

if mol.tags['step'] < 2:

    #Single points with TDDFT
    os.makedirs('SP',exist_ok=True)
    os.chdir('SP')

    mol.RefineConformers(GAUSSIAN,
                         jobname='{0}-SP'.format(input_name),
                         runtype='sp',
                         method=dftmethod,
                         scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                         nproc=nproc,
                         mem=mem)

    os.chdir('../')

    mol.tags['step'] = 2
    mol.checkpoint('{0}.chk'.format(input_name))
    mol.checkpoint('{0}-SP.mol'.format(input_name))

##################################################################

if mol.tags['step'] < 3:

    mol.conformers = mol.conformers[0:10]

    oldchk = None
    geom = None
    guess = None

    #check if the chk exists from a previous run
    if os.path.exists('lowest10/{0}-opt-0/{0}-opt-0-try0.chk'.format(input_name)):
        prefix = '{0}-opt-'.format(input_name)
        oldchk = prefix + '{conf}/previous.chk'
        geom = 'check'
        guess = 'read'

        #check which try was the last one for each
        #rename the last run file to previous.chk
        for index in range(0,10):
            if os.path.exists('lowest10/{0}-opt-{1}/{0}-opt-{1}-try3.chk'.format(input_name,index)):
                os.rename('lowest10/{0}-opt-{1}/{0}-opt-{1}-try3.chk'.format(input_name,index),'lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            elif os.path.exists('lowest10/{0}-opt-{1}/{0}-opt-{1}-try2.chk'.format(input_name,index)):
                os.rename('lowest10/{0}-opt-{1}/{0}-opt-{1}-try2.chk'.format(input_name,index),'lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            elif os.path.exists('lowest10/{0}-opt-{1}/{0}-opt-{1}-try1.chk'.format(input_name,index)):
                os.rename('lowest10/{0}-opt-{1}/{0}-opt-{1}-try1.chk'.format(input_name,index),'lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            elif os.path.exists('lowest10/{0}-opt-{1}/{0}-opt-{1}-try0.chk'.format(input_name,index)):
                os.rename('lowest10/{0}-opt-{1}/{0}-opt-{1}-try0.chk'.format(input_name,index),'lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            else:
                raise IndexError('Could not find chk to resubmit')


    #Optimize lowest10 conformers
    os.makedirs('lowest10',exist_ok=True)
    os.chdir('lowest10')

    mol.RefineConformers(GAUSSIAN,
                         jobname='{0}-opt'.format(input_name),
                         runtype='opt_freq',
                         opt='ts,noeigen,recalcfc=10',
                         method=dftmethod,
                         scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                         nproc=nproc,
                         mem=mem,
                         oldchk=oldchk,
                         geom=geom,
                         guess=guess,
                         tries=3)

    os.chdir('../')

    #write out the TS energy from the lowest conf
    df = pd.DataFrame([[mol.conformers[0].properties['free_energy']]],columns=['TS'])
    df.to_csv(output_energies,index=False)

    mol.tags['step'] = 3
    mol.checkpoint('{0}.chk'.format(input_name))
    mol.checkpoint('{0}-lowest10.mol'.format(input_name))

##################################################################

if mol.tags['step'] < 4:

    #take the lowest energy conformer
    mol = mol.conformers[0]

    mol.ToXYZ('{0}-lowest_conf.xyz'.format(input_name))

    os.makedirs('IRC',exist_ok=True)
    os.chdir('IRC')

    calc = GAUSSIAN(mol,
                   runtype='irc',
                   jobname='{0}-IRC-forward'.format(input_name),
                   method=dftmethod,
                   scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                   irc='forward,calcfc',
                   nproc=nproc,
                   mem=mem)

    forward = Run(calc,tries=1)

    os.chdir('../')

    df = pd.read_csv(output_energies)
    df['forward'] = [forward.energy]
    df.to_csv(output_energies,index=False)

    forward.ToXYZ('{0}-forward.xyz'.format(input_name))

    mol.tags['step'] = 4
    mol.checkpoint('{0}.chk'.format(input_name))
    forward.checkpoint('{0}-forward.mol'.format(input_name))

#######################################################################

if mol.tags['step'] < 5:

    os.chdir('IRC')

    calc = GAUSSIAN(mol,
                   runtype='irc',
                   jobname='{0}-IRC-reverse'.format(input_name),
                   method=dftmethod,
                   scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                   irc='reverse,calcfc',
                   nproc=nproc,
                   mem=mem)

    reverse = Run(calc,tries=1)

    os.chdir('../')

    df = pd.read_csv(output_energies)
    df['reverse'] = [reverse.energy]
    df.to_csv(output_energies,index=False)

    reverse.ToXYZ('{0}-reverse.xyz'.format(input_name))

    mol.tags['step'] = 5
    mol.checkpoint('{0}.chk'.format(input_name))
    reverse.checkpoint('{0}-reverse.mol'.format(input_name))

#######################################################################

if mol.tags['step'] < 6:

    if not forward:
        forward = molecule.ReadCheckpoint('{0}-forward.mol'.format(input_name))
    if not reverse:
        reverse = molecule.ReadCheckpoint('{0}-reverse.mol'.format(input_name))

    #determine which is the product by comparing atom distances
    def distance(x1,y1,z1,x2,y2,z2):
        return(np.sqrt(   ((x1-x2)**2 ) + ((y1-y2)**2 ) + ((z1-z2)**2 ) ))

    #check the constrained distances to see which is larger
    forward_coords1 = forward.coords[atom1]
    forward_coords2 = forward.coords[atom2]
    forward_distance = distance(forward_coords1[0],
                                forward_coords1[1],
                                forward_coords1[2],
                                forward_coords2[0],
                                forward_coords2[1],
                                forward_coords2[2])

    reverse_coords1 = reverse.coords[atom3]
    reverse_coords2 = reverse.coords[atom4]
    reverse_distance = distance(reverse_coords1[0],
                                reverse_coords1[1],
                                reverse_coords1[2],
                                reverse_coords2[0],
                                reverse_coords2[1],
                                reverse_coords2[2])

    df = pd.read_csv(output_energies)
    df['forward_distance'] = forward_distance
    df['reverse_distance'] = reverse_distance

    if forward_distance > reverse_distance:
        product = reverse
        df['product_direction'] = ['reverse']
    else:
        product = forward
        df['product_direction'] = ['forward']

    df.to_csv(output_energies,index=False)

    mol.tags['step'] = 6
    mol.checkpoint('{0}.chk'.format(input_name))
    product.checkpoint('{0}-product.mol'.format(input_name))

###############################################################

if mol.tags['step'] < 7:

    if not product:
        product = molecule.ReadCheckpoint('{0}-product.mol'.format(input_name))

    oldchk = None
    geom = None
    guess = None


    #check if the chk exists from a previous run
    if os.path.exists('product/{0}-product/{0}-product-try0.chk'.format(input_name)):
        oldchk = 'product/{0}-product/previous.chk'.format(input_name)
        geom = 'check'
        guess = 'read'

        #check which try was the last one for each
        #rename the last run file to previous.chk
        if os.path.exists('product/{0}-product/{0}-product-try3.chk'.format(input_name)):
            os.rename('product/{0}-product/{0}-product-try3.chk'.format(input_name),'product/{0}-product/previous.chk'.format(input_name))
        elif os.path.exists('product/{0}-product/{0}-product-try2.chk'.format(input_name)):
            os.rename('product/{0}-product/{0}-product-try2.chk'.format(input_name),'product/{0}-product/previous.chk'.format(input_name))
        elif os.path.exists('product/{0}-product/{0}-product-try1.chk'.format(input_name)):
            os.rename('product/{0}-product/{0}-product-try1.chk'.format(input_name),'product/{0}-product/previous.chk'.format(input_name))
        elif os.path.exists('product/{0}-product/{0}-product-try0.chk'.format(input_name)):
            os.rename('product/{0}-product/{0}-product-try0.chk'.format(input_name),'product/{0}-product/previous.chk'.format(input_name))
        else:
            raise IndexError('Could not find chk to resubmit')


    os.makedirs('product',exist_ok=True)
    os.chdir('product')

    calc = GAUSSIAN(reverse,
                   runtype='opt_freq',
                   jobname='{0}-product'.format(input_name),
                   method=dftmethod,
                   scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                   nproc=nproc,
                   oldchk=oldchk,
                   geom=geom,
                   check=check,
                   mem=mem)

    product = Run(calc,tries=3)

    os.chdir('../')

    df = pd.read_csv(output_energies)
    df['product'] = [product.energy]
    df.to_csv(output_energies,index=False)

    product.ToXYZ('{0}-product.xyz'.format(input_name))

    mol.tags['step'] = 7
    mol.checkpoint('{0}.chk'.format(input_name))
    product.checkpoint('{0}-product-optimized.mol'.format(input_name))

################################################################

if mol.tags['step'] < 8:

    if not product:
        product = molecule.ReadCheckpoint('{0}-product-optimized.mol'.format(input_name))

    os.makedirs('product-conf_search',exist_ok=True)
    os.chdir('product-conf_search')



    conformer_calculator = CREST(product,
                                 jobname='{0}-confsearch'.format(input_name),
                                 runtype='confsearch',
                                 nproc=nproc,
                                 mem=mem,
                                 gbsa=solvent_dict['crest'],
                                 max_confs=100)

    mol = Run(conformer_calculator)

    os.chdir('../')

    mol.tags['step'] = 8
    mol.checkpoint('{0}.chk'.format(input_name))
    product.checkpoint('{0}-product-confsearch.mol'.format(input_name))

##################################################################

if mol.tags['step'] < 9:

    if not product:
        product = molecule.ReadCheckpoint('{0}-product-confsearch.mol'.format(input_name))

    #Single points with TDDFT
    os.makedirs('product-SP',exist_ok=True)
    os.chdir('product-SP')

    product.RefineConformers(GAUSSIAN,
                         jobname='{0}-SP'.format(input_name),
                         runtype='sp',
                         method=dftmethod,
                         scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                         nproc=nproc,
                         mem=mem)

    os.chdir('../')

    mol.tags['step'] = 9
    mol.checkpoint('{0}.chk'.format(input_name))
    product.checkpoint('{0}-product-SP.mol'.format(input_name))

##################################################################

if mol.tags['step'] < 10:

    if not product:
        product = molecule.ReadCheckpoint('{0}-product-SP.mol'.format(input_name))

    product.conformers = product.conformers[0:10]

    oldchk = None
    geom = None
    guess = None

    #check if the chk exists from a previous run
    if os.path.exists('product-lowest10/{0}-opt-0/{0}-opt-0-try0.chk'.format(input_name)):
        prefix = '{0}-opt-'.format(input_name)
        oldchk = prefix + '{conf}/previous.chk'
        geom = 'check'
        guess = 'read'

        #check which try was the last one for each
        #rename the last run file to previous.chk
        for index in range(0,10):
            if os.path.exists('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try3.chk'.format(input_name,index)):
                os.rename('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try3.chk'.format(input_name,index),'product-lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            elif os.path.exists('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try2.chk'.format(input_name,index)):
                os.rename('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try2.chk'.format(input_name,index),'product-lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            elif os.path.exists('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try1.chk'.format(input_name,index)):
                os.rename('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try1.chk'.format(input_name,index),'product-lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            elif os.path.exists('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try0.chk'.format(input_name,index)):
                os.rename('product-lowest10/{0}-opt-{1}/{0}-opt-{1}-try0.chk'.format(input_name,index),'product-lowest10/{0}-opt-{1}/previous.chk'.format(input_name,index))
            else:
                raise IndexError('Could not find chk to resubmit')


    #Optimize lowest10 conformers
    os.makedirs('product-lowest10',exist_ok=True)
    os.chdir('product-lowest10')

    product.RefineConformers(GAUSSIAN,
                         jobname='{0}-opt'.format(input_name),
                         runtype='opt_freq',
                         method=dftmethod,
                         scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                         nproc=nproc,
                         mem=mem,
                         oldchk=oldchk,
                         geom=geom,
                         guess=guess,
                         tries=3)

    os.chdir('../')

    #write out the TS energy from the lowest conf
    df = pd.read_csv(output_energies)
    df['global_minima_product'] = [product.conformers[0].properties['free_energy']]
    df.to_csv(output_energies,index=False)

    mol.tags['step'] = 10
    mol.checkpoint('{0}.chk'.format(input_name))
    product.checkpoint('{0}-product-lowest10.mol'.format(input_name))

    product.conformers[0].ToXYZ('{0}-product-global-min.xyz'.format(input_name))


##################################################################

print('Done!')

