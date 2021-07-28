from workflow import molecule
from workflow.crest import CREST
from workflow.gaussian import GAUSSIAN
from workflow.calculator import Run,RunBatch

from workflow.utils import parallelize
import os


def get_lowest_conformer(smiles,name=None):
     ########################################################################
    # STEP 1 #Initialize a molecule object for the smiles NCC(CC1CCCC1O)=O #
    ########################################################################
    if name is None:
        #check if tuple of smiles,name is passed
        try:
            smiles, name = smiles
        except:
            raise TypeError("No 'name' input found: either explicitly pass the 'name' argument, or pass a tuple of smiles,name")

    mol = molecule.SmilesToMol(smiles,nconfs=1)

    ########################################################
    # STEP 2 #Generate CREST conformers for the mol object #
    ########################################################
    #mol.constraints = [(1,2),(5,9),(1,2,3),(4,5,6,7)]

    os.makedirs('conf_search',exist_ok=True)
    os.chdir('conf_search')

    crest_conformer_calculator = CREST(mol,title='{0}-crest'.format(name),runtype='confsearch',nproc=16,mem=120,time='1-00:00:00',partition='short,lopez',gbsa='water',ewin=500)

    mol = Run(crest_conformer_calculator)

    os.chdir('../')

    ##################################################################
    # STEP 3 #Re-rank the conformers by Gaussian single point energy #
    ##################################################################

    os.makedirs('SP',exist_ok=True)
    os.chdir('SP')

    mol.refine_conformers(GAUSSIAN,nproc=16,mem=120,jobname='{0}-sp-refine-conformers'.format(name),runtype='sp',method='B3LYP/6-31G',empiricaldispersion='GD3bj')

    os.chdir('../')

    ###################################################################
    # STEP 4 #Optimize the lowest 10 conformers and rank by Free energy
    ###################################################################

    os.makedirs('OPT',exist_ok=True)
    os.chdir('OPT')

    #if there are less then 10, optimize them all
    if len(mol.conformers) > 10:
        mol.conformers = mol.conformers[0:10]

    mol.refine_conformers(GAUSSIAN,nproc=16,mem=120,jobname='{0}-opt-refine-conformers'.format(name),runtype='opt_freq',method='B3LYP/6-31G',empiricaldispersion='GD3bj',opt='calcfc')

    os.chdir('../')

    #################################################
    # STEP 5 #Write out the lowest energy conformer #
    #################################################

    os.makedirs('lowest_energy_conformers',exist_ok=True)
    os.chdir('lowest_energy_conformers')

    mol.conformers[0].ToXYZ('{0}-lowestconf.xyz'.format(name))

    os.chdir('../')

    return((mol,name))

def get_lowest_conformers(smiles_name_pairs,nproc):
    results = parallelize(smiles_name_pairs,get_lowest_conformer,nproc=nproc)
    return(results)
