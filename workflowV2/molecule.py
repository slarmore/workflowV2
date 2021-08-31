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
import pickle

#######################################################################################################
#class for freezing the attributes of mol and conformer objects
class FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and key not in dir(self):
            raise TypeError( "Cannot add attribute to {0}, use .tags dict to store misc. information".format(self))
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True



#######################################################################################################
#Mol object

class Mol(FrozenClass):
    '''Main molecule class'''
    def __init__(self,
                atoms,
                coords,
                smiles=None,
                conformers=[],
                energy=None,
                constraints=[],
                tags={},
                charge=0,
                mult=1,
                properties={}, 
                warnings=[]): 


    #######################
    #check for valid input#
        if len(atoms) != len(coords):
            raise IndexError('Length of atoms is {0}, length of coords is {1}'.format(len(atoms),len(coords)))

    #################
    #core attributes#
        self._coords = False
        self._atoms = atoms
        self.coords = coords
        self.conformers = conformers
        self.energy = energy
        self.constraints = constraints
        self.tags = tags
        self.warnings = warnings
        self.charge = charge
        self.mult = mult
        self.properties = properties
        self._natoms = len(self.atoms)
        self.smiles = smiles


    ####################
    #derived attributes#
        self.get_rdkitmol()


        #don't allow any new attribute to be defined
        self._freeze()

    def get_xyz(self):
        if len(self.atoms) != len(self.coords):
            raise IndexError('Length of atoms is {0}, length of coords is {1}'.format(len(self.atoms),len(self.coords)))

        xyz = []

        for index,atom in enumerate(self.atoms):
            #need to prevent writing out scientific notation by specifying number of decimals
            xyz.append([atom,'{:.7f}'.format(self.coords[index][0]),'{:.7f}'.format(self.coords[index][1]),'{:.7f}'.format(self.coords[index][2])])
        self.xyz = np.array(xyz)

    def get_rdkitmol(self):
        if self.smiles is not None:
            try:
                self.rdkitmol = Chem.MolFromSmiles(self.smiles)
            except:
                warning('Could not generate RDkit mol from {0}'.format(self.smiles))
        else:
            self.rdkitmol = None
        
    def get_xyzstring(self):
        xyzstring=[]
        for line in self.xyz:
            xyzstring.append('      '.join(line))
        self.xyzstring = '\n'.join(xyzstring)

    def update_geometry(self):
        self.get_xyz()
        self.get_xyzstring()

    #restrict access to some Mol properties that user
    #shouldn't touch, or ensure updates to the rest of the info
    @property
    def atoms(self):
        return(self._atoms)

    @atoms.setter
    def atoms(self,new_atoms):
        self._atoms = np.array(new_atoms)
        self._natoms = len(new_atoms)
        if self._coords:
            self.update_geometry()
            self._coords = True

    @property
    def coords(self):
        return(self._coords)

    @coords.setter
    def coords(self,new_coords):
        if len(new_coords) != len(self.atoms):
            raise IndexError('Length of atoms is {0}, length of coords is {1}'.format(len(self.atoms),len(self.coords)))

        #make a numpy array of float points for consistency
        self._coords = np.array(new_coords,dtype=np.float64)
        self.update_geometry()

    @property
    def natoms(self):
        return(self._natoms)

###################
#File IO functions#
    def ToXYZ(self,outfile=None,title=None,write=True):
        #number of atoms
        out = [str(self.natoms)]

        #title line
        if title is not None:
            out.append(str(title))
        elif 'title' in self.tags:
            out.append(str(self.tags['title']))
        else:
            out.append(str(self.energy))

        out.append(self.xyzstring)
        out.append('\n')

        #write to file if requrested, otherwise return string
        if write:
            if outfile is None:
                raise TypeError('outfile must not be None if write is True')
                
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

                out.append('\n')

            with open(outfile,'w') as xyzout:
                xyzout.write('\n'.join(out))
            
        else:
            raise IndexError('No conformers present')

    def checkpoint(self,outfile):
        with open(outfile,'wb') as chk:
            pickle.dump(self,chk)

###########################################################################
#copy method to prevent changing original object with updated calculations#
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
                charge=deepcopy(self.charge),
                mult=deepcopy(self.mult),
                properties=deepcopy(self.properties), 
                warnings=deepcopy(self.warnings))
        return(_copy)

    def copy(self):
        return(deepcopy(self))

#######################
#convinience functions#
    def RefineConformers(self,calculator_class,jobname,tries=1,ignore=False,**kwargs):
        '''take in a generic calculator to apply to each conformer
        then re-rank the conformers with relative energies'''

        #take the generic calculator and apply it to each conformre
        #calculators = [calculator_class(conf,jobname=jobname+'-{0}'.format(index),**kwargs) for index,conf in enumerate(self.conformers)]
        
        #look for {conf} in any of the kwarg values and replace with conformer #
        #this will allow specification of specific oldchk file for each conformer
        #when restarting a gaussian calculation
        calculators = []
        for index,conf in enumerate(self.conformers):
            local_kwargs = kwargs.copy()
            for key, value in local_kwargs.items():
                if  isinstance(value,str):
                    local_kwargs[key] = value.replace('{conf}',str(index))
            calculators.append(calculator_class(conf,jobname=jobname+'-{0}'.format(index),**local_kwargs))


        conformers = calculator.RunBatch(calculators,jobname,tries=tries,ignore=ignore)
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



#######################################################################################################
#Conformer object

class Conformer(FrozenClass):
    '''a container to store conformers - its a lighter weight version of the mol object'''
    def __init__(self,
                    atoms,
                    coords,
                    energy=None,
                    tags={},
                    constraints=[],
                    smiles=None,
                    charge=0,
                    mult=1):

        #######################
        #check for valid input#
        if len(atoms) != len(coords):
            raise IndexError('Length of atoms is {0}, length of coords is {1}'.format(len(atoms),len(coords)))

        self._coords = False
        self._atoms = atoms
        self.coords = coords
        self.energy = energy
        self.constraints = constraints
        self.tags = tags
        self.charge = charge
        self.mult = mult
        self._natoms = len(self.atoms)
        self.smiles = smiles

        #don't allow any new attribute to be defined
        self._freeze()

    def get_xyz(self):
        if len(self.atoms) != len(self.coords):
            raise IndexError('Length of atoms is {0}, length of coords is {1}'.format(len(self.atoms),len(self.coords)))

        xyz = []

        for index,atom in enumerate(self.atoms):
            #need to prevent writing out scientific notation by specifying number of decimals
            xyz.append([atom,'{:.7f}'.format(self.coords[index][0]),'{:.7f}'.format(self.coords[index][1]),'{:.7f}'.format(self.coords[index][2])])
        self.xyz = np.array(xyz)
        
    def get_xyzstring(self):
        xyzstring=[]
        for line in self.xyz:
            xyzstring.append('      '.join(line))
        self.xyzstring = '\n'.join(xyzstring)

    def update_geometry(self):
        self.get_xyz()
        self.get_xyzstring()

    #restrict access to some Conformer properties that user
    #shouldn't touch, or ensure updates to the rest of the info
    @property
    def atoms(self):
        return(self._atoms)

    @atoms.setter
    def atoms(self,new_atoms):
        self._atoms = new_atoms
        self._natoms = len(new_atoms)
        self.update_geometry()

    @property
    def coord(self):
        return(self.coord)

    @coord.setter
    def coord(self,new_coords):
        if len(new_coords) != len(self.atoms):
            raise IndexError('Length of atoms is {0}, length of coords is {1}'.format(len(self.atoms),len(self.coords)))

        #make a numpy array of float points for consistency
        self.coords = np.array(new_coords,dtype=np.float64)
        log(type(self.coords[0][0]))
        self.update_geometry()

    @property
    def natoms(self):
        return(self._natoms)

###################
#File IO functions#
    def ToXYZ(self,outfile,title=None,write=True):
        #number of atoms
        out = [str(self.natoms)]

        #title line
        if title is not None:
            out.append(str(title))
        elif 'title' in self.tags:
            out.append(str(self.tags['title {0}'.format(self.energy)]))
        else:
            out.append(str(self.energy))

        out.append(self.xyzstring)

        #write to file if requrested, otherwise return string
        if write:
            with open(outfile,'w') as xyzout:
                xyzout.write('\n'.join(out))
        else:
            return('\n'.join(out))

#######################
#upgrade to mol object#
    def ToMol(self):
        mol = Mol(self.atoms,
                self.coords,
                smiles=self.smiles,
                conformers=[],
                energy=self.energy,
                constraints=self.constraints,
                tags=self.tags,
                charge=self.charge,
                mult=self.mult,
                properties={})
        return(mol)



###########################################################################################################################
# Type converters

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
                    constraints=constraints,tags=tags,
                    charge=charge,mult=mult))
    else:
        warning('No such file {0}'.format(infile))
        return(None)

def SmilesToMol(smiles,
            nconfs=5,
            tags={},
            charge=0,
            mult=1,
            constraints=[],
                
            #aruguments for embedmultipleconfs
            randomSeed=0,
            maxAttempts=5000,
            numThreads=0,
            pruneRmsThresh=0.005,
            useRandomCoords=False,
            enforceChirality=True,
            boxSizeMult=2.0,
            randNegEig=True,
            numZeroFail=1,
            coordMap={},
            orceTol=0.001,
            ignoreSmoothingFailures=False,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
            printExpTorsionAngles=False,
            useSmallRingTorsions=False,
            useMacrocycleTorsions=False,
         
            #arguments for conformer opts
            opt='UFF',
            maxIters=10000,
            mmffVarient='MMFF94',
            nonBondedThresh=100.0,
            ignoreInterfragInteractions=True,
            vdwThresh=10.0):
    '''Wrapper for taking in a smiles and using RDKIT to generate conformers and energies'''

    if nconfs < 1:
        raise IndexError('Must generate at least 1 conformer')
    if not opt in ['UFF','MMFF']:
        raise KeyError("opt argument must be either 'UFF' or 'MMF'")

    rdkitmol = Chem.MolFromSmiles(smiles)
    rdkitmol = AllChem.AddHs(rdkitmol,addCoords=True)
    
    num_generated = 0
    while num_generated < nconfs:
        confs = AllChem.EmbedMultipleConfs(rdkitmol,
                               randomSeed=randomSeed,
                               maxAttempts=maxAttempts,
                               numThreads=numThreads,
                               nonBondedThresh=nonBondedThresh,
                               ignoreInterfragInteractions=ignoreInterfragInteractions,
                               pruneRmsThresh=pruneRmsThresh,
                               useRandomCoords=useRandomCoords,
                               enforceChirality=enforceChirality,
                               boxSizeMult=boxSizeMult,
                               randNegEig=randNegEig,
                               numZeroFail=numZeroFail,
                               coordMap=coordMap,
                               forceTol=forceTol,
                               ignoreSmoothingFailures=ignoreSmoothingFailures,
                               useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
                               useBasicKnowledge=useBasicKnowledge,
                               printExpTorsionAngles=printExpTorsionAngles,
                               useSmallRingTorsions=useSmallRingTorsions,
                               useMacrocycleTorsions=useMacrocycleTorsions)
        
        num_generated = len(confs)
        pruneRmsThresh -= -0.001
        log('lowered rms threshold to {0} to find more than {1} conformers'.format(pruneRmsThresh,num_generated))
    
    if opt == 'UFF':
        energies = AllChem.UFFOptimizeMoleculeConfs(rdkitmol,
                                                    maxIters=maxIters,
                                                    ignoreInterfragInteractions=ignoreInterfragInteractions,
                                                    vdwThresh=vdwThresh)
    elif opt == 'MMFF':
        energies = AllChem.MMFFOptimizeMoleculeConfs(rdkitmol,
                                                     maxIters=maxIters,
                                                     mmffVarient=mmffVarient,
                                                     nonBondedThresh=nonBondedThresh,
                                                     ignoreInterfragInteractions=ignoreInterfragInteractions)
    else:
        raise KeyError("opt argument must be either 'UFF' or 'MMF'")
            
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
        conformers.append(Conformer(energy=energy,atoms=np.array(atoms),coords=np.array(coords),tags=tags,smiles=smiles,charge=charge,mult=mult))

    #sort the conformers by energy (first value in the tuple)
    conformers.sort(key=lambda x: x.energy)

    #use the lowest conformer energy as the definition for the mol object
    coords = conformers[0].coords
    atoms = conformers[0].atoms
    energy = conformers[0].energy

    #return the mol object
    return(Mol(atoms=atoms,coords=coords,smiles=smiles,energy=energy,conformers=conformers,
                    constraints=constraints,tags=tags,
                    charge=charge,mult=mult))


def ReadCheckpoint(infile):
    with open(infile,'rb') as chk:
        return(pickle.load(chk))

