from workflow import molecule
from workflow.gaussian import GAUSSIAN
from workflow.calculator import Run

mol = molecule.XYZToMol('badmol.xyz')

calculator = GAUSSIAN(mol,runtype='opt_freq',nproc=16,mem=120,jobname='opt_badmol-{0}'.format(tries),opt='calcfc,ts,noeigen,maxstep=10',method='b3lyp/6-31G(d)',empiricaldispersion='GD3BJ',scrf='iefpcm,solvent=water',scf='xqc')

mol = Run(calculator)

print('Optimized the geometry successfully')

print('test')
print('test again')