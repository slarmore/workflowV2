#!/usr/bin/env python
#SBATCH --job-name=TS_workflow
#SBATCH --output=out-%A
#SBATCH --error=error-%A
#SBATCH --partition=lopez,long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=4-00:00:00

solvent_dict = {'crest':'DMSO','gaussian':'dimethylsulfoxide'}
dftmethod = 'M062X/6-31G'
nproc = 12
mem = 60
time = '1-00:00:00'

###################################################################

from workflowV2 import molecule
from workflowV2.utils import atom_distance
from workflowV2.calculator import Run as Run
from workflowV2.software.GAUSSIAN import GAUSSIAN as GAUSSIAN
from workflowV2.software.CREST import CREST
from workflowV2 import message
import pandas as pd
import sys
import os
import numpy as np
from datetime import datetime
from argparse import ArgumentParser

###################################################################

def setup(input_file,
         charge,
         multiplicity,
         bond1,
         bond2):

    '''Get things prepared before jumping into workflow'''
    global atom1,atom2,atom3,atom4,input_name,output_energies,mol,forward,reverse,product


    atom1 = bond1[0]
    atom2 = bond1[1]
    atom3 = bond2[0]
    atom4 = bond2[1]

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
        message.log('''##################################################################
    found checkpoined mol for restart - {0}
    ##################################################################

    '''.format(datetime.now()),time=False)
        message.log('restarting at step {0}'.format(mol.tags['step']+1))


    if not mol:
        message.log('''##################################################################

    Starting workflow at - {0}

    ##################################################################

    '''.format(datetime.now()),time=False)
        mol = molecule.XYZToMol(input_file,charge=charge,mult=multiplicity)
        mol.constraints = [[atom1,atom2],[atom3,atom4]]

        mol.tags['step'] = 0

##################################################################

def main(input_file,
         charge,
         multiplicity,
         bond1,
         bond2):

    '''Multi-step TS workflow'''

    #Get some global variables set up before starting the flow
    #not the nicest way to do so ... I assume no one is importing this main function outside of the script
    setup(input_file,charge,multiplicity,bond1,bond2)

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


        os.chdir('../')

        mol.tags['step'] = 1
        mol.checkpoint('{0}.chk'.format(input_name))
        mol.checkpoint('{0}-confsearch.mol'.format(input_name))

    ##################################################################

    if mol.tags['step'] < 2:

        #Single points to re-rank conformers
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

        #Optimize lowest10 conformers
        os.makedirs('lowest10',exist_ok=True)
        os.chdir('lowest10')

        mol.conformers = mol.conformers[0:10]

        mol.RefineConformers(GAUSSIAN,
                            jobname='{0}-preopt'.format(input_name),
                            runtype='opt',
                            method=dftmethod,
                            scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                            nproc=nproc,
                            mem=mem,
                            tries=1,
                            time=time)

        mol.tags['step'] = 3
        mol.checkpoint('{0}.chk'.format(input_name))
        mol.checkpoint('{0}-preopt.mol'.format(input_name))

        os.chdir('../')

        mol.tags['step'] = 3
        mol.checkpoint('{0}.chk'.format(input_name))
        mol.checkpoint('{0}-preopt.mol'.format(input_name))

    ##################################################################

    if mol.tags['step'] < 4:

        #release the constraints now that I have the conformers I needed
        mol.constraints = []
        for conformer in mol.conformers:
            conformer.constraints = []

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
                            tries=3,
                            TS=True,
                            time=time)

        os.chdir('../')

        mol.tags['step'] = 4
        mol.checkpoint('{0}.chk'.format(input_name))
        mol.checkpoint('{0}-lowest10.mol'.format(input_name))

    ##################################################################

    if mol.tags['step'] < 5:


        #check that the structure still resembles a TS
        #we know the bond lengths for the forming bonds should be less than 3A
        #take the first one that is a TS on to the next step

        valid_TS = False
        for index,conformer in enumerate(mol.conformers):
            distance1 = atom_distance(conformer,atom1,conformer,atom2)
            distance2 = atom_distance(conformer,atom3,conformer,atom4)

            #if distances are both between 2 and 3A, then it's likely a TS
            if 2.0 < distance1 < 3.0  and 2.0 < distance2 < 3.0 :
                message.log('Conformer{0} has TS bond distances {1}, {2} - VALID TS found'.format(index,distance1,distance2))
                valid_TS = True
                mol.tags['valid_TS'] = index
                break
            else:
                message.log('Conformer{0} has TS bond distances {1}, {2} - NOT valid TS'.format(index,distance1,distance2))

        if not valid_TS:
            raise IndexError('No valid TSs')

        #write out a csv of the conformer energies
        conformerdf = pd.DataFrame([['conformer{0}'.format(index),conformer.properties['free_energy'],conformer.energy] for index,conformer in enumerate(mol.conformers)],columns=['Conformer','Free_energy','Relative_energy'])
        conformerdf.to_csv('{0}-TS_conformers.csv'.format(input_name),index=False)

        #write out the TS energy from the lowest valid TS
        df = pd.DataFrame([[mol.conformers[valid_TS].properties['free_energy']]],columns=['TS_free_energy'])
        df.to_csv(output_energies,index=False)

        mol.tags['step'] = 5
        mol.checkpoint('{0}.chk'.format(input_name))
        mol.checkpoint('{0}-lowest10.mol'.format(input_name))


    ##################################################################

    if mol.tags['step'] < 6:

        #take the lowest energy valid TS
        #read from the tags in case this is a restarted run
        mol = mol.conformers[mol.tags['valid_TS']]

        mol.ToXYZ('{0}-lowest_conf.xyz'.format(input_name))

        os.makedirs('IRC',exist_ok=True)
        os.chdir('IRC')

        calc = GAUSSIAN(mol,
                    runtype='irc',
                    jobname='{0}-IRC-forward'.format(input_name),
                    method=dftmethod,
                    scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                    irc='forward,rcfc',
                    oldchk=mol.tags['chk'],
                    nproc=nproc,
                    mem=mem,
                    time=time)

        forward = Run(calc,tries=1,ignore=True)

        os.chdir('../')

        #check if at least 3 steps were taken
        if len(forward.properties['optimization_energies']) < 3:
            message.warning('IRC didnt take any steps')
            exit()

        df = pd.read_csv(output_energies)
        df['IRC_forward_electronic_energy'] = [forward.energy]
        df.to_csv(output_energies,index=False)

        forward.ToXYZ('{0}-forward.xyz'.format(input_name))

        mol.tags['step'] = 6
        mol.checkpoint('{0}.chk'.format(input_name))
        forward.checkpoint('{0}-forward.mol'.format(input_name))

    #######################################################################

    if mol.tags['step'] < 7:

        os.chdir('IRC')

        calc = GAUSSIAN(mol,
                    runtype='irc',
                    jobname='{0}-IRC-reverse'.format(input_name),
                    method=dftmethod,
                    scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                    irc='reverse,rcfc',
                    oldchk=mol.tags['chk'],
                    nproc=nproc,
                    mem=mem,
                    time=time)

        reverse = Run(calc,tries=1,ignore=True)

        os.chdir('../')

        #check if at least 3 steps were taken
        if len(reverse.properties['optimization_energies']) < 3:
            message.warning('IRC didnt take any steps')
            exit()


        df = pd.read_csv(output_energies)
        df['IRC_reverse_electronic_energy'] = [reverse.energy]
        df.to_csv(output_energies,index=False)

        reverse.ToXYZ('{0}-reverse.xyz'.format(input_name))

        mol.tags['step'] = 7
        mol.checkpoint('{0}.chk'.format(input_name))
        reverse.checkpoint('{0}-reverse.mol'.format(input_name))

    #######################################################################

    if mol.tags['step'] < 8:

        if not forward:
            forward = molecule.ReadCheckpoint('{0}-forward.mol'.format(input_name))
        if not reverse:
            reverse = molecule.ReadCheckpoint('{0}-reverse.mol'.format(input_name))

        #check the constrained distances to see which is larger
        forward_distance = atom_distance(forward,atom1,forward,atom2)
        reverse_distance = atom_distance(reverse,atom1,reverse,atom2)

        df = pd.read_csv(output_energies)
        df['forward_distance'] = forward_distance
        df['reverse_distance'] = reverse_distance

        if forward_distance > reverse_distance:
            product = reverse.copy()
            df['product_direction'] = ['reverse']
        else:
            product = forward.copy()
            df['product_direction'] = ['forward']

        df.to_csv(output_energies,index=False)

        mol.tags['step'] = 8
        mol.checkpoint('{0}.chk'.format(input_name))
        product.checkpoint('{0}-product.mol'.format(input_name))

    ###############################################################

    if mol.tags['step'] < 9:

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

        calc = GAUSSIAN(product,
                    runtype='opt_freq',
                    jobname='{0}-product'.format(input_name),
                    method=dftmethod,
                    scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                    nproc=nproc,
                    oldchk=oldchk,
                    geom=geom,
                    guess=guess,
                    mem=mem,
                    time=time)

        product = Run(calc,tries=3)

        os.chdir('../')

        df = pd.read_csv(output_energies)
        df['product_local_minima_free_energy'] = [product.energy]
        df.to_csv(output_energies,index=False)

        product.ToXYZ('{0}-product.xyz'.format(input_name))

        mol.tags['step'] = 9
        mol.checkpoint('{0}.chk'.format(input_name))
        product.checkpoint('{0}-product-optimized.mol'.format(input_name))

    ################################################################

    if mol.tags['step'] < 10:

        if not product:
            product = molecule.ReadCheckpoint('{0}-product-optimized.mol'.format(input_name))

        os.makedirs('product-conf_search',exist_ok=True)
        os.chdir('product-conf_search')



        conformer_calculator = CREST(product,
                                    jobname='{0}-product-confsearch'.format(input_name),
                                    runtype='confsearch',
                                    nproc=nproc,
                                    mem=mem,
                                    gbsa=solvent_dict['crest'],
                                    max_confs=100,
                                    time=time)

        product = Run(conformer_calculator)

        os.chdir('../')

        mol.tags['step'] = 10
        mol.checkpoint('{0}.chk'.format(input_name))
        product.checkpoint('{0}-product-confsearch.mol'.format(input_name))

    ##################################################################

    if mol.tags['step'] < 11:

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
                            mem=mem,
                            time=time)

        os.chdir('../')

        mol.tags['step'] = 11
        mol.checkpoint('{0}.chk'.format(input_name))
        product.checkpoint('{0}-product-SP.mol'.format(input_name))

    ##################################################################

    if mol.tags['step'] < 12:

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
                            jobname='{0}-product-opt'.format(input_name),
                            runtype='opt_freq',
                            method=dftmethod,
                            scrf='iefpcm,solvent={0}'.format(solvent_dict['gaussian']),
                            nproc=nproc,
                            mem=mem,
                            oldchk=oldchk,
                            geom=geom,
                            guess=guess,
                            tries=3,
                            time=time)

        os.chdir('../')

        #write out a csv of the conformer energies
        conformerdf = pd.DataFrame([['conformer{0}'.format(index),conformer.properties['free_energy'],conformer.energy] for index,conformer in enumerate(product.conformers)],columns=['Conformer','Free_energy','Relative_energy'])
        conformerdf.to_csv('{0}-product_conformers.csv'.format(input_name),index=False)

        #write out the product energy from the lowest conf
        df = pd.read_csv(output_energies)
        df['global_minima_product_free_energy'] = [product.conformers[0].properties['free_energy']]
        df.to_csv(output_energies,index=False)

        mol.tags['step'] = 12
        mol.checkpoint('{0}.chk'.format(input_name))
        product.checkpoint('{0}-product-lowest10.mol'.format(input_name))

        product.conformers[0].ToXYZ('{0}-product-global-min.xyz'.format(input_name))


if __name__ == '__main__':

    #parse arguments
    parser = ArgumentParser('TS workflow with 2 constrained bonds')
    parser.add_argument('-i','--input',type=str,dest='input_file')
    parser.add_argument('-q','--charge',type=int,dest='charge',default=0)
    parser.add_argument('-m','--multiplicity',type=int,dest='multiplicity',default=1)
    parser.add_argument('--bond1',nargs=2,type=int)
    parser.add_argument('--bond1',nargs=2,type=int)
    arguments = parser.parse_args()

    main(arguments.input_file,
         arguments.charge,
         arguments.multiplicity,
         arguments.bond1,
         arguments.bond2)

    message.log('Done!')
