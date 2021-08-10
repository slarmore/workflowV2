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
 And is created with the bash command:\
 `conda env create -f workflowV2.yml`

Alternatively, you can install directly with pip, **HOWEVER** RDKIT and numpy are dependencies that the user must manage themselves:\
`pip install git+https://github.com/neal-p/workflowV2.git`

To use any calculators and actually compute anything, you will need to provide the paths to the executables in the `config.py` file and may need to change the default slurm parameters. 

For exmple, for Discovery users at Northeastern the G16 root is set to `g16root = '/work/lopez/'`, which points to our installation of Gaussian16.\
Similarly, to run CREST, `crest_exe = '/work/lopez/xtb/crest'` and `xtb_exe = '/work/lopez/xtb/xtb_6.2.3/bin/xtb'` are specified. 

Additional global configuration variables are also accessible here. For example, the default Slurm partition is set `default_partition = 'short,lopez'`

# Updating
Updates can be pulled with the following:\
` pip  install --upgrade git+https://github.com/neal-p/workflowV2.git`

# Creating and using Mol Objects

# Creating and using Calculators
