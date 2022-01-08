# Conformational Search Routine

Performs CREST conformational search, followed by SinglePoint energy re-ranking, and full optimization of X lowest energy conformers

## Output
Lowest energy conformer xyz

Lowest X conformers combined xyz

csv containing CREST, SinglePoint, fully optimized energies, and final relative energies


## In this directory
Script to run is called `conf_search_routine.sbatch`, and takes in an xyz file of the molecule to conformational search, as well as a number of arguments to control the workflow

A test molecule file `testmol.xyz` is provided, as well as the file containing the arguments used for an example run (`args.txt`)

The results of an example run are in the `SampleResults` folder. 

## Running

#### Pre-requisits
  - an xyz file of the molecule to conformational search
  - up-to-date workflowV2_env environment active
 
 #### Execution
 The script can be called on the command line and given the input xyz by `sbatch conf_search_routine --xyz testmol.xyz`
 
 
 However, more often then not, you'll want to modify the method, solvent, resources, etc. The following arguments are availible:
 
 ```
 usage: Generalized Conformational Search Routine [-h] [--xyz MOL] [-o LOGFILE] [-n NAME] [--maxconfs MAX_CONFS] [-q CHARGE] [-m MULT] [-b BONDS BONDS] [-a ANGLES ANGLES ANGLES] [-d DIHEDRALS DIHEDRALS DIHEDRALS DIHEDRALS] [--TS] [--final_opts N_FINAL_OPTS] [--solvent SOLVENT]
                                                 [--CREST_solvent CREST_SOLVENT] [--method METHOD] [--nproc NPROC] [--mem MEM] [--partition PARTITION] [--time TIME]

optional arguments:
  -h, --help            show this help message and exit
  --xyz MOL             molecule input file in xyz format
  -o LOGFILE, --logfile LOGFILE
                        name of file to log workflow output, defaults to guessing a name based on the xyz file
  -n NAME, --name NAME  molecule name for calculation input/output files etc, defaults to guessing a name based on the xyz file
  --maxconfs MAX_CONFS  maximim number of conformers to retrieve from CREST
  -q CHARGE, --charge CHARGE
                        molecule charge
  -m MULT, --mult MULT  molecule multiplicity
  -b BONDS BONDS, --bond BONDS BONDS
                        2 atom indicies to constrain a distance, 0 indexed - only applied to CREST conformational search, not subsequent optimizations
  -a ANGLES ANGLES ANGLES, --angle ANGLES ANGLES ANGLES
                        3 atom indicies to constain an angle, 0 indexed - only applied to CREST conformational search, not subsequent optimizations
  -d DIHEDRALS DIHEDRALS DIHEDRALS DIHEDRALS, --dihedral DIHEDRALS DIHEDRALS DIHEDRALS DIHEDRALS
                        4 atom indicies to constrain a dihedral, 0 indexed - only applied to CREST conformational search, not subsequent optimizations
  --TS                  Flag to treat the input as a transition state structure
  --final_opts N_FINAL_OPTS
                        Maximum number of conformers to keep for the final full optimization
  --solvent SOLVENT     Solvent to use for calculations, if not specified will run in gas phase. See https://gaussian.com/scrf/ and https://xtb-docs.readthedocs.io/en/latest/gbsa.html?highlight=solvent for availible solvents
  --CREST_solvent CREST_SOLVENT
                        Separate solvent parameter to use for CREST conformational search only - can be set to "None". To be used if the final desired solvent is not availible in CREST, for example. See https://xtb-docs.readthedocs.io/en/latest/gbsa.html?highlight=solvent for availible solvents.
  --method METHOD       Calculation method (dft functional and basis set for example)
  --nproc NPROC         Number of processors to use for calculations
  --mem MEM             Memory to be allocated to each calculation (Gb)
  --partition PARTITION
                        Partition(s) to run calculations on
  --time TIME           Time limit for each individual calculation
```

Since specifying a number of these arguments on the comnmand-line can be cumberson at best, and woefully error prone at worst, I prefer to list my arguments out in a file. You can then tell Python to read arguments from the file by using the argument file prefix `@`:


args.txt >>>
```
--xyz=testmol.xyz
--dihedral
10
9
8
2
--method=B3LYP/6-31G(d,p)
--solvent=water
--CREST_solvent=None
--nproc=16
--mem=120
--final_opts=10
```

Submitted with
`sbatch conf_search_routine.sbatch @args.txt`

This example will run a CREST conformational search on the molecule in testmol.xyz in gas phase with the 10,9,8,2 atom dihedral constrained, then re-rank the conformers by B3LYP/6-31G(d,p) IEFPCM(water) SinglePoint energy, and finally fully optimize the 10 lowest energy conformers with B3LYP/6-31G(d,p) IEFPCM(water)

Top of the output csv is shown below:
![Example_csv](example_csv.png)



