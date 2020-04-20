from ase.calculators.turbomole import Turbomole
from ase.io import read

### normal mode calculation with optimized structure, just zncluster

mol = read('znclusir_pbe.xyz')

params = {'task': 'normal mode analysis',
    'density convergence': 1.0e-7}

calc = Turbomole(restart = True,**params)

mol.set_calculator(calc)
calc.calculate(mol)

results = calc.get_results()


f = open('results.txt', 'r') 

print(calc.runtime, f)
print(calc['results']['vibrational spectrum'],f)
print(calc.todict(skip_default = False),f)

f.close()






