# GAUSSIAN calculation interface
Submits gaussian calculations on the input molecule

The gaussian calculator object is created with the GAUSSIAN function

```
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2 import molecule
from workflowV2.calculator import Run

mol = molecule.SmilesToMol('C1=CC=CC=C1')
calc = GAUSSIAN(mol,runtype='opt_freq',method='HF/6-31G') 
result = Run(calc)
```

## Gaussian calculation options

### Required arguments
`mol` - the molecule object to compute\
`jobname` - name of the calculation (this names the directory and files for the calculation)\
`runtype` - what type of calculation is to be run:
  - `'sp'`
  - `'opt'`
  - `'freq'`
  - `'opt_freq'`
  - `'irc'`

`method` - the gaussian method to use\
  for example, a dft calculation should give the `functional/basis set`\
  a semi-empirical calculation should give the name of the semi-empirical method like `'PM7'`

### Job related arguments
`nproc` - number of processors\
`mem` - memory in GB\
`time` - time limit for the job\
`partition` - partition for the job

### Calculation specifications
Most gaussian keywords for the route line are availible for use [gaussian manual](http://gaussian.com/man/)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **General formatting** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Route line options can be specified like:\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `calc = GAUSSIAN(mol,runtype='opt_freq',method='B3LYP/6-31G',opt='calcfc')`\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; where the `keyword = 'some option'` is added to the route line as `keyword=(some option)`
&nbsp; \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Another example would be to add implicit solvent with gaussian's scrf keyword\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `calc = GAUSSIAN(mol,runtype='opt_freq',method='HF/6-31G',opt='calcfc',scrf='iefpcm,solvent=toluene')`\
&nbsp; \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Finally, we can add empirical dispersion correction with
```
calc = GAUSSIAN(mol,runtype='opt_freq',opt='calcfc'
                   method='HF/6-31G',empiricaldispersion='GD3bj',scrf='iefpcm,solvent=toluene'
                   )
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (it's often useful to split up the function call to multiple lines since complex calcuations can have lengthy route lines)\
&nbsp; \
&nbsp; \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Dealing with Constraints** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If the input Mol object has constraints defined in the `.constraints` attribute, \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the GAUSSIAN calculator will automatically add `modredundant` to the `opt` keyword,  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; apply the appropriate atom constraints after the input section. \
&nbsp; \
&nbsp; \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Reading from chk** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Additionally, keywords for reading from chk files are allowed, but the `oldchk` should be specified \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `calc = GAUSSIAN(mol,runtype='opt_freq',method='HF/6-31G',oldchk='/home/benzene.chk',geom='check',guess='read')` \
&nbsp; \
&nbsp; \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Advanced options** \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Lastly, there are two special keywords for advanced use: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - `afterinput` - This allows the user to specify additional input for after the input coordinates. \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This can be used for specifying custom solvents \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - `advanced_route_opts` - This is a catch-all keyword that will put anything passed to it directly on the route line, exactly as written: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `advanced_route_opts='banana formcheck'` is added to the route line as `banana formcheck` 

