#This script is the base class for the mol object

from rdkit import Chem
from rdkit.Chem import AllChem
import os
import numpy as np
from copy import copy, deepcopy
from . import calculator 
from .utils import hartree2kcal
import itertools
from .message import warning,log,display
from .config import max_conformers

class Mol:
    '''Main molecule class'''
    def __init__(self,
                atoms,
                coords,
                smiles=None,
                conformers=[],
                energy=None,
                constraints=[],
                tags={},
                xyz=None,
                charge=0,
                mult=1,
                properties={}, 
                calc=None,
                warnings=[]): 

        #core properties
        self.smiles = smiles
        self.coords = coords
        self.conformers = conformers
        self.energy = energy
        self.constraints = constraints
        self.tags = tags
        self.warnings = warnings
        self.atoms = atoms
        self.charge = charge
        self.mult = mult
        self.properties = properties
        self.calc = calc
        self.natoms = len(self.atoms)

        #get the derived properties
        if xyz is None:
            self.get_xyz()
        else:
            self.xyz = xyz

        self.get_rdkitmol()
        self.get_xyzstring()
        

    #derived property methods
    def get_xyz(self):
        xyz = []
        for index,atom in enumerate(self.atoms):
            xyz.append([atom,self.coords[index][0],self.coords[index][1],self.coords[index][2]])
        self.xyz = np.array(xyz)

    def get_rdkitmol(self):
        if self.smiles is not None:
            self.rdkitmol = Chem.MolFromSmiles(self.smiles)
        else:
            self.rdkitmol = None
        
    def get_xyzstring(self):
        xyzstring=[]
        for line in self.xyz:
            #need to prevent writing out scientific notation by specifying number of decimals
            line = ['{:.7f}'.format(value) if type(value) == np.float64 else value for value in line]
            xyzstring.append('      '.join(line))
        self.xyzstring = '\n'.join(xyzstring)

    #Mol object methods

    def update_geometry(self):
        self.get_xyz()
        self.get_xyzstring()


    #File IO functions
    def ToXYZ(self,outfile=None,title=None,coords=None,write=True):
        #number of atoms
        out = [str(self.natoms)]

        #title line
        if title is not None:
            out.append(str(title))
        elif 'title' in self.tags:
            out.append(str(self.tags['title {0}'.format(self.energy)]))
        else:
            out.append(str(self.energy))

        #pair the atom with the coordinates
        if coords is None:
            coords = self.coords
        for atom,line in zip(self.atoms,list(coords)):
            line = [str(coord) for coord in line]
            line.insert(0,atom)
            out.append('      '.join(line))

        #write to file if requrested, otherwise return string
        if write:
            if outfile is None:
                raise TypeError('outfile argument cannot be None')

            with open(outfile,'w') as xyzout:
                xyzout.write('\n'.join(out))
        else:
            return('\n'.join(out))

    def ConformersToXYZ(self,outfile,start=0,stop=None):
        if len(self.conformers) > 0:
            if stop is None:
                stop = len(self.conformers) + 1
            out = []
            for i,conf in enumerate(self.conformers[start:stop]):
                if i > 0:
                    out.append('')
                out_sub = [str(conf.natoms)]
                if 'title' in self.tags:
                    out_sub.append(str('{0} - {1}'.format(self.tags['title'],conf.energy)))
                else:
                    out_sub.append('From Conformer object - {0}'.format(conf.energy))
                for line in conf.xyz:
                    out_sub.append('      '.join(line))
                out.append('\n'.join(out_sub))

            with open(outfile,'w') as xyzout:
                xyzout.write('\n'.join(out))
            
        else:
            raise IndexError('No conformers present')

    
    def __deepcopy__(self,memo):
        id_self = id(self)        # memoization avoids unnecesary recursion
        _copy = memo.get(id_self)
        if _copy is None:
            _copy = type(self)(deepcopy(self.atoms),
                deepcopy(self.coords),
                smiles=deepcopy(self.smiles),
                conformers=deepcopy(self.conformers),
                energy=deepcopy(self.energy),
                constraints=deepcopy(self.constraints),
                tags=deepcopy(self.tags),
                xyz=deepcopy(self.xyz),
                charge=deepcopy(self.charge),
                mult=deepcopy(self.mult),
                properties=deepcopy(self.properties), 
                calc=deepcopy(self.calc),
                warnings=deepcopy(self.warnings))
        return(_copy)

    def copy(self):
        return(deepcopy(self))

    def RefineConformers(self,calculator_class,jobname,**kwargs):
        '''take in a generic calculator to apply to each conformer
        then re-rank the conformers with relative energies'''

        #take the generic calculator and apply it to each conformre
        calculators = [calculator_class(conf,jobname=jobname+'{0}'.format(index+1),**kwargs) for index,conf in enumerate(self.conformers)]
        
        conformers = calculator.RunBatch(calculators,jobname)
        conformers.sort(key=lambda x: x.energy)
        energies = np.array([conf.energy for conf in conformers])
        energies = energies - energies[0]
        energies = np.vectorize(hartree2kcal)(energies)
        
        #assign back on the mol
        self.conformers = conformers
        for conf,energy in zip(self.conformers,energies):
            conf.energy = energy

    def generate_rdkit_conformers(self,nconfs=5,tags={}):
        if self.rdkitmol is None:
            if self.smiles is None:
                raise TypeError('No SMILES for this mol object')
            else:
                self.rdkitmol = Chem.MolFromSmiles(self.smiles)
        self.rdkitmol = AllChem.AddHs(self.rdkitmol,addCoords=True)
        AllChem.EmbedMultipleConfs(self.rdkitmol,numConfs=nconfs,enforceChirality=True,numThreads=0)
        energies = AllChem.UFFOptimizeMoleculeConfs(self.rdkitmol,numThreads=0)
        energies = [energy[1] for energy in energies]


        if 'origin' not in tags:
            tags['origin'] = 'rdkit'

        #store a conformer object for each conformer
        #that should make it easy to turn them into full Mol objects later if needed
        conformers =[]
        for conf in range(0,nconfs):
            coords = []
            atoms = []
            energy = energies[conf]
            for atom in range(0,self.rdkitmol.GetNumAtoms()):
                symbol = self.rdkitmol.GetAtomWithIdx(atom).GetSymbol()
                x = self.rdkitmol.GetConformer(conf).GetAtomPosition(atom).x
                y = self.rdkitmol.GetConformer(conf).GetAtomPosition(atom).y
                z = self.rdkitmol.GetConformer(conf).GetAtomPosition(atom).z
                
                coords.append([x,y,z])
                atoms.append(symbol)
            conformers.append(Conformer(energy,np.array(atoms),np.array(coords),xyz=None,tags=tags,smiles=self.smiles,charge=self.charge,mult=self.mult))

        #sort the conformers by energy (first value in the tuple)
        conformers.sort(key=lambda x: x.energy)


class Conformer:
    '''a container to store conformers - its a lighter weight version of the mol object'''
    def __init__(self,energy,atoms,coords,xyz=None,tags={},constraints=[],smiles=None,charge=0,mult=1):
        self.energy = energy
        self.atoms = atoms
        self.coords = coords
        self.tags = tags
        self.natoms = len(atoms)
        self.constraints = constraints
        self.smiles = smiles
        self.charge = charge
        self.mult = mult

        #construct the xyz if not present
        if xyz is None:
            xyz = []
            for index,atom in enumerate(atoms):
                xyz.append([atom,coords[index][0],coords[index][1],coords[index][2]])
            self.xyz = np.array(xyz)
        else:
            self.xyz = xyz

        xyzstring=[]
        for line in self.xyz:
            xyzstring.append('      '.join(line))
        self.xyzstring = '\n'.join(xyzstring)
    
        #File IO functions
    def ToXYZ(self,outfile,title=None):
        out = [str(self.natoms)]
        if title is not None:
            out.append(str(title))
        elif 'title' in self.tags:
            out.append(str(self.tags['title']))
        else:
            out.append('From Mol object - {0}'.format(self.energy))
        for line in self.xyz:
            out.append('      '.join(line))

        with open(outfile,'w') as xyzout:
            xyzout.write('\n'.join(out))
    
    def ToMol(self):
        mol = Mol(self.atoms,
                self.coords,
                smiles=self.smiles,
                conformers=[],
                energy=self.energy,
                constraints=self.constraints,
                tags=self.tags,
                xyz=self.xyz,
                charge=self.charge,
                mult=self.mult,
                properties={}, 
                calc=None)
        return(mol)

###########################################################################################################################
# Type converters
##########################################################################################################################

def XYZToMol(infile,
            smiles=None,
            energy=None,
            constraints=[],
            tags={},
            charge=0,
            mult=1):

    if os.path.exists(infile):
        with open(infile,'r') as xyzfile:
            rawxyz = xyzfile.read().splitlines()
            
            #first line should be number of atoms
            try:
                natom = int(rawxyz[0])
            except:
                #need to learn how to raise these things better
                warning('First line of {0} is not an integer\n It should be the number of atoms'.format(infile))
                return(None)
            title = rawxyz[1]
            rawxyz = rawxyz[2:2+natom]

        coords = []
        atoms = []
        xyz = []

        for line in rawxyz:
            atom,x,y,z = line.split()
            xyz.append([atom,x,y,z])
            atoms.append(atom)
            coords.append([x,y,z])
        
        xyz = np.array(xyz)
        atoms = np.array(atoms)
        coords = np.array(coords)

        if not 'title' in tags:
            tags['title'] = title 


        #return the mol object
        return(Mol(atoms,coords,smiles=smiles,energy=energy,
                    constraints=[],tags=tags,
                    charge=charge,mult=mult,
                    xyz=xyz))
    else:
        warning('No such file {0}'.format(infile))
        return(None)

def SmilesToMol(smiles,
            nconfs=5,
            tags={},
            charge=0,
            mult=1):
    '''Wrapper for taking in a smiles and using RDKIT to generate conformers and energies'''
   
    if nconfs > max_conformers:
        raise ValueError('The max_conformers is set to {0}'.format(max_conformers))

    rdkitmol = Chem.MolFromSmiles(smiles)
    rdkitmol = AllChem.AddHs(rdkitmol,addCoords=True)
    AllChem.EmbedMultipleConfs(rdkitmol,numConfs=nconfs,enforceChirality=True,numThreads=0)
    energies = AllChem.UFFOptimizeMoleculeConfs(rdkitmol,numThreads=0)
    energies = [energy[1] for energy in energies]

    if 'origin' not in tags:
        tags['origin'] = 'rdkit'
    
    #store a conformer object for each conformer
    #that should make it easy to turn them into full Mol objects later if needed
    conformers =[]
    for conf in range(0,nconfs):
        coords = []
        atoms = []
        energy = energies[conf]
        for atom in range(0,rdkitmol.GetNumAtoms()):
            symbol = rdkitmol.GetAtomWithIdx(atom).GetSymbol()
            x = rdkitmol.GetConformer(conf).GetAtomPosition(atom).x
            y = rdkitmol.GetConformer(conf).GetAtomPosition(atom).y
            z = rdkitmol.GetConformer(conf).GetAtomPosition(atom).z
            
            coords.append([x,y,z])
            atoms.append(symbol)
        conformers.append(Conformer(energy,np.array(atoms),np.array(coords),xyz=None,tags=tags,smiles=smiles,charge=charge,mult=mult))

    #sort the conformers by energy (first value in the tuple)
    conformers.sort(key=lambda x: x.energy)

    #use the lowest conformer energy as the definition for the mol object
    coords = conformers[0].coords
    atoms = conformers[0].atoms
    energy = conformers[0].energy

    #return the mol object
    return(Mol(atoms,coords,smiles=smiles,energy=energy,conformers=conformers,
                    constraints=[],tags=tags,
                    charge=charge,mult=mult,
                    xyz=None))





