# workflow2.0 Overview

A package for automating QM calculation workflows

Inspired by RDkit and the Atomic Simulation Environment's design principles\
but adapted for HPC job scheduler environment

More here soon!



## The Mol Object
The central piece to this package is the **Mol Object**. 

<img src="Mol_Object.png" width="75%" align="middle" />

The Mol Object stores all of the 3D information, propreties, conformers, energies, tags, etc for a molecule. The central theme of this codebase is manipulating these Mol Objects with calculators that will modify the Mol Object in some way. The Mol Object can be passed to a calculator, which will return an updated Mol Object. 

## Calculators
Submit some kind of calculation to update the Mol Object

<img src="Calculator.png" width="75%" align="middle"/>

For example, creating a Gaussian calculator to optimize the input Mol Object would return a Mol Object with updated coordinates, energies, orbital energies, vibrational modes, etc. 

## Creating Workflows
This code is designed to build quantum chemistry workflows.


<table>
<col style="width:50%">
<col style="width:50%">
<thead>
<tr>
<td style="text-align: left; vertical-align:top;"> <p> The power of this codebase is not for individual  <br />calculations.  Often those are easy enough to set up <br /> and manage by hand.  Rather, it's utility is in the <br />ability to script together complex workflows, passing the Mol Object form one calculator to the next, <br />automating many calculations. </p> By building out functionality from the central <br />Mol Object, everything will connect nicely together <br />and allow much more flexibility when designing automated pipelines compared to stringing together sbatch scripts and hoping for the best. </td>
<td> <img src="Workflow.png" width="150%"/> </td>
</tr>
</thead>
</table>

# Installation
The easiest way to install this code is by creating a new conda environment with the following yml file:
``` 
name: workflowV2_env
channels:
  - conda-forge
  - defaults
dependencies:
  - rdkit
  - pip
  - pip:
    - git+https://github.com/neal-p/workflowV2.git
 ```
 Create the environment with the bash command:\
 `conda env create -f workflowV2.yml`

Alternatively, you can install directly with pip, **HOWEVER** RDKIT and numpy are dependencies that the user then must manage themselves:\
`pip install git+https://github.com/neal-p/workflowV2.git`

To use any calculators and actually compute anything, you will need to provide the paths to the software executables in the `config.py` file and may need to change the default slurm parameters. 

For exmple, for Discovery users at Northeastern the G16 root is set to `g16root = '/work/lopez/'`, which points to our installation of Gaussian16.\
Similarly, to run CREST, `crest_exe = '/work/lopez/xtb/crest'` and `xtb_exe = '/work/lopez/xtb/xtb_6.2.3/bin/xtb'` are specified. 

Additional global configuration variables are also accessible here. For example, the default Slurm partition is set `default_partition = 'short,lopez'`

# Updating
Updates can be pulled with the following:\
` pip  install --upgrade git+https://github.com/neal-p/workflowV2.git`

# Creating and using Mol Objects
The `molecule` module has all of tools to create and manage Mol Objects

The Mol Object typically isn't created directly, so this documentation will focus on the attributes and methods for using the Mol Objects

```
from workflowV2 import molecule
mol = molecule.SmilesToMol('C1=CC=CC=C1')
```

  - `mol.atoms` - np.array of atomic symbols
  - `mol.coords` - np.array of x y z coordinates
  - `mol.conformers` - List of Conformer Objects (a lighter-weight Mol object)
  - `mol.energy` - molecule energy
  - `mol.constraints` - list of lists containing constraints: [[1],[1,2],[1,2,3],[1,2,3,4]] will constain atom1, the distance between atom1-atom2, the angle atom1-atom2-atom3, and the dihedral atom1-atom2-atom3-atom4 
  - `mol.tags` - A dictionary of tags for a molecule. This can be helpful for checkpointing long workflows with `mol.tags['step'] = 3` or adding information like `mol.tags['global_min'] = True`
  - `mol.warnings` - 
  - `mol.charge` - 
  - `mol.mult` - 
  - `mol.properties` - 
  - `mol.natoms` - 
  - `mol.smiles` - 
  - `mol.rdkitmol` - 
  - `mol.xyz` - 
  - `mol.xyzstring` - 
  - `mol.ToXYZ(outfile=None,title=None,write=True)` - 
  - `mol.ConformersToXYZ(outfile,start=0,stop=None)` - 
  - `mol.checkpoint(outfile)` - 
  - `mol.copy()` - 
  - `mol.RefineConformers(calculator_class,jobname,tries=1,**kwargs)` - 
  - `mol.generate_rdkit_conformers(nconfs=5,tags={})` - 




## From SMILES
Mol Objects can be created from SMILES strings with the `SmilesToMol` function. 

`molecule.SmilesToMol(smiles,nconfs=5,tags={},charge=0,mult=1,constraints=[],seed=0)`
  - `smiles` - input SMILES string
  - `nconfs` - number of conformers to generate, minimum value is 1
  - `tags` - initialize dictionary of tags for the mol
  - `charge` - molecule charge
  - `mult` - molecule multiplicity
  - `constraints` - initialize constraints for the mol
  - `seed` - random seed for 3D conformer generation by RDkit

Initializing Mol Object:
```
from workflowV2 import molecule
mol = molecule.SmilesToMol('C1=CC=CC=C1')
```

## From XYZ file
Mol Objects can be created from XYZ file with the `XYZToMol` function. 

`molecule.XYZToMol(infile,smiles=None,energy=None,constraints=[],tags={},charge=0,mult=1)`
  - `infile` - input XYZ file
  - `smiles` - smiles string
  - `energy` - molecule energy
  - `constraints` - initialize constraints for the mol
  - `tags` - initialize dictionary of tags for the mol
  - `charge` - molecule charge
  - `mult` - molecule multiplicity

Example XYZ file:
```
USER@example:/home$ cat benzene.xyz

>6
>benzene coordinates
>
>
>
>
>
>
>
>
>
>
>
>

```

Initializing Mol Object:
```
from workflowV2 import molecule
mol = molecule.XYZToMol('benzene.xyz')
```

# Creating and using Calculators
