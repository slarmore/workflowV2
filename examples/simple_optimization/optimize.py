from workflowV2 import molecule
from workflowV2.calculator import Run as Run
from workflowV2.software.GAUSSIAN import GAUSSIAN 
from workflowV2 import message
from argparse import ArgumentParser


def main(smiles,charge=0,multiplicity=1):
    '''This program will create 3D coordinates for a SMILES string
     and preform a very simple HF/6-31G optimization with Gaussian.
     
     The user can optionally set the charge and multiplicity for the calculation.
     
     The optimized XYZ coordinates and energy will be written to 
         - optimized.xyz
         - energy.txt
    '''
    
  
    #Create Mol Object
    mol = molecule.SmilesToMol(smiles,charge=charge,mult=multiplicity)
    
    #Create Gaussian calculator
    calc = GAUSSIAN(mol,jobname='simple_optimization-job',runtype='opt',method='HF/6-31G')
    
    #Run calculator
    result = Run(calc)
    
    #write optimized xyz coordinates
    result.ToXYZ('optimized.xyz')
    
    #write the energy
    with open('energy.txt','w') as energyfile:
        energyfile.write('The electronic energy is {0}'.format(result.energy))
        
    
if __name__ == '__main__':
  
    #parse smiles, charge, multiplicity
    parser = ArgumentParser('Optimize with HF/6-31G')
    parser.add_argument('-s','--smiles',type=str,dest='smiles')
    parser.add_argument('-q','--charge',type=int,dest='charge',default=0)
    parser.add_argument('-m','--multiplicity',type=int,dest='multiplicity',default=1)
    arguments = parser.parse_args()
  
    #optimize the molecule
    print('\n############################\nSetting up for optimization!\n\n...\n\n')
    main(arguments.smiles,charge=arguments.charge,multiplicity=arguments.multiplicity)
  
    #Let the user know it's done
    message.log('\nDone optimizing!')
