from ase import Atoms
import numpy as np
from ase.calculators.turbomole import Turbomole
from ase.constraints import FixBondLength
import csv
'''
    Calculation of energies and forces of 20 different water-configuration
'''
# Turbomole params set up
params = {
        'title': 'water',
        'use redundant internals': True,
        'basis set name': 'def2-SV(P)',
        'total charge': 0,
        'multiplicity': 1,
        'use dft': True,
        'density functional': 'b3-lyp',
        'use resolution of identity': True,
         }
# Create 20 different values for H-O Distance and H-O-H angle
r0 = 0.969
roH = np.arange(0.9 ,2.0, 0.055, dtype = "float")
theta  = np.arange(0.5,1.0, 0.025)
# Dicts to save energy and forces for every water configuration
energies = {}
forces = {}
# Creation of 20 different water molecule geometries
for r in roH:
    for t in theta:
    
        pos = [(r0, 0, 0),  # H1
               (r * np.cos(t), r * np.sin(t),0), # H2
               (0, 0, 0)] # O
        
        atoms = Atoms('H2O', positions=pos)
        atoms.constraints = FixBondLength(0,2)
    
        a = 8.0 # cell parameter
        h = 0.2 # grid space
    
        atoms.set_cell((a,a,a))
        atoms.center()       
    # Turbomole calc set up
        atoms.calc=Turbomole(**params)
        energy = atoms.get_potential_energy()
        force = atoms.get_forces()
        dictkey = str(r)+"_"+str(t)
        energies.update({dictkey:energy}) 
        forces.update({dictkey:force})
        
# Write System energies and forces to csv file
w = csv.writer(open("energy_output.csv","w"))
for key, val in energies.items():
	w.writerow([key,val])
w.close()
v = csv.writer(open("force_output.csv","w"))
for key, val in forces.items():
	w.writerow([key,val])
v.close()