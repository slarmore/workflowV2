from workflow.workflows import get_lowest_conformers

compounds = [('O=C(C(CCO)N)C1=CC=CC=C1','mol1'),('CCCCCCCCCCCCCCC','mol2'),('O=C1C(CC2=CC=CC(Br)=C2C)CCCC1','mol3'),('NC(C1=CC=CC=C1)(O)CC2CC2','mol3')]

#calls python multiprocessing to run in paraellel
results = get_lowest_conformers(compounds,nproc=4)

energies = [result.energy for result in results]

print(energies)
