# workflow2.0

A package for automating QM calculation workflows\

Inspired by RDkit and the Atomic Simulation Environment's design principles\
but adapted for HPC job scheduler environment\

More here soon!



# The Mol Object
The central piece to this package is the **Mol Object**. 
![alt text](Mol_Object.png)

The Mol Object stores all of the 3D information, propreties, conformers, energies, tags, etc for a molecule. The central theme of this codebase is manipulating these Mol Objects with calculators that will modify the Mol Object in some way. The Mol Object can be passed to a calculator, which will return an updated Mol Object. 

![alt text](Calculator.png)

For example, creating a Gaussian calculator to optimize the input Mol Object would return a Mol Object with updated coordinates, energies, orbital energies, vibrational modes, etc. 
