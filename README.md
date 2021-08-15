# workflow2.0 Overview

A package for automating QM calculation workflows

Inspired by RDkit and the Atomic Simulation Environment's design principles\
but adapted for a QM focus in the HPC job scheduler environment

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

**Basic properties and methods**
  - `mol.atoms` (np.array) -  atomic symbols
  - `mol.coords` (np.array) - x y z coordinates
  - `mol.charge` (int) - moleucle charge
  - `mol.mult` (int) - molecule multiplicity
  - `mol.natoms` (int) - number of atoms
  - `mol.smiles` (string) - SMILES string (if provided, or created from smiles)
  - `mol.rdkitmol` (object) - RDkit molecule if SMILES is present
  - `mol.energy` (float) - molecule energy
  - `mol.tags` (dict) - tags for a molecule. This can be helpful for checkpointing long workflows with `mol.tags['step'] = 3` or adding information like `mol.tags['global_min'] = True`
  - `mol.checkpoint(outfile)` - pickles the Mol Object and writes to `outfile`
  - `mol.constraints` (list of lists) - apply atom constraints for any calculations. Each sub-list can contain 1 atom (will constrain xyz), two atoms (constrain the distance), three atoms (constrain the angle), or four atoms (constarin the dihedral): `[[1],[1,2],[1,2,3],[1,2,3,4]]` will constain atom1, the distance between atom1-atom2, the angle atom1-atom2-atom3, and the dihedral atom1-atom2-atom3-atom4 
  - `mol.copy()` (Mol Object) - return a deep copy of the Mol Object 

**Related to coordinates and coordinate files**
  - `mol.xyz` (np.array) - xyz coordinates
  - `mol.xyzstring` (string) - string representation of the xyz coordinates
  - `mol.ToXYZ(outfile=None,title=None,write=True)` - write coordinates to an xyz file called `outfile`, with `title` line, unless `write` is False, in which case a string of the xyz file is returned

**Related to conformers**
  - `mol.conformers` (list) - Conformer Objects (a lighter-weight Mol object)
  - `mol.generate_rdkit_conformers(nconfs=5,tags={})` - Uses RDkit's `embed_multiple_confs` to generate conformers
  - `mol.ConformersToXYZ(outfile,start=0,stop=None)` - writes a multi-structure xyz file to `output` with the Mol Object's conformers from index `start` to `stop`
  - `mol.RefineConformers(calculator_class,jobname,tries=1,**kwargs)` - create a calculator, runs it for every conformer associated with the Mol Object, and re-rank the conformers by the updated energies. The original Mol Object is updated, there is no return value. The type of calculation is specified with the `calculator_class`, a general `jobname` that will be used for the calculations, and the number of `tries` should be specified to determine if auto-resubmissions will be attempted. All other arguments are exactly the same as running the calculation on a single Mol Object. For example, a typical Gaussian calculation is submitted like 
 ```
 calc = GAUSSIAN(mol,jobname='benzene_opt',runtype='opt_freq',method='HF/6-31G')
 ```
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; to apply the same calculation to each conformer would look like
 ```
 mol.RefineConformers(GAUSSIAN,jobname='benzene_opt',tries=3,runtype='opt_freq',method='HF/6-31G')
 ```

**Related to calculation results**
  - `mol.warnings` (list) - Warnings generated by calculations. For exmaple, the `SCF_Error` warning will be added if the SCF doesn't converge in a Gaussian calculation
  - `mol.properties` (dict) - Properties parsed after a calculation will be added to the properties dictionary. For example, after a Gaussian frequency calculation, the free energy correction is stored at `mol.properties['free_energy_correction']`. See the specific software documentation for the key of your desired property



## Creating from SMILES
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

## Creating from XYZ file
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
>C      -1.2250455      0.6633388      -0.1287105
>C      -1.1950034      -0.7271697      0.0226333
>C      0.0300421      -1.3905085      0.1513437
>C      1.2250455      -0.6633388      0.1287106
>C      1.1950034      0.7271697      -0.0226332
>C      -0.0300421      1.3905085      -0.1513439
>H      -2.1727752      1.1765163      -0.2282845
>H      -2.1194917      -1.2897286      0.0401430
>H      0.0532835      -2.4662450      0.2684274
>H      2.1727752      -1.1765163      0.2282849
>H      2.1194917      1.2897287      -0.0401430
>H      -0.0532835      2.4662449      -0.2684279

```

Initializing Mol Object:
```
from workflowV2 import molecule
mol = molecule.XYZToMol('benzene.xyz')
```

# Creating and using Calculators
The `calculator` module has the functions for interacting with Calculator Objects that are created by the program-specific files in `workflowV2.software.XXX`

The process of running calculations has two steps: creating the calculator and submitting the calculator

**Creating the calculator** \
A Calculator Object is created to run a specific program. Each supported program has a file in the workflowV2/software/ directory and the function to create a Calculator Object is all-caps the name of the program. For example, to set up a Gaussian calculator you need to: 
  1. import the GAUSSIAN function from workflowV2.software.GAUSSIAN
  2. call the GAUSSIAN function with the proper arguments to set up the desired calculation

```
from workflowV2.software.GAUSSIAN import GAUSSIAN
calc = GAUSSIAN(mol,jobname='benzene_opt',runtype='opt_freq',method='HF/6-31G')
```
**Submitting the calculator** \
The submission process is handled by the funcitons in the `calculator` module. There are two ways of submitting:
  - `calculator.Run` - submit a single calculator 
  - `calculator.Runbatch` - submits a list of calculators in parallel with a Slurm array

Both methods return updated Mol Objects from the calculation. The `Run` funtion returns a single Mol Object back, the `RunBatch` function returns the list of updated Mol Objects.

`calculator.Run` arguments:
  - `calculator` - the calculator object to submit the job with
  - `tries=1` - the number of resubmissions if the job fails (this is only supported for programs where the `fix_errors` method is implemented)
  - `ignore=False` - ignore calculation errors on the final resubmission try. When False, if the calculation fails after all of the tries are used up, it will return the updated Mol Object with whatever information could be parsed. For example, even if a structure did not finish optimizing, it will return the latest geometry. When False, if the final try fails an IndexError (`Ran out of resubmission tries`) will be raised.


`calculator.RunBatch` arguments:
  - `calculators` - a list of calculator objects to submit the job with
  - `jobname='batch_job'` - Name of the Slurm array job to submit
  - `max=50` - number of concurrent array jobs to run
  - `tries=1` - the number of resubmissions if the job fails (this is only supported for programs where the `fix_errors` method is implemented)
  - `ignore=False` - ignore calculation errors on the final resubmission try. When False, if the calculation fails after all of the tries are used up, it will return the updated Mol Object with whatever information could be parsed. For example, even if a structure did not finish optimizing, it will return the latest geometry. When False, if the final try fails an IndexError (`Ran out of resubmission tries`) will be raised.

# utils functions
I have built a few utility functions that I've found helpful along the way in `workflowV2.utils` that include things for parallel computing, atom distances, etc. Hopefully these will be encorporated more smoothly into the rest of the code as it evolves!

# Example workflows
There are a few workflow exmaples of varrying complexity in the `examples` directory

# Reporting bugs and requests
Please report any bugs or feature requests! This hopefully a growing and expanding code base that will adapt to our needs!
