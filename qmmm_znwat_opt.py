#Created on Tue Nov 19 09:39:39 2019

#@author: tombolelli



import numpy as np
from ase.calculators.tip3p import TIP3P, epsilon0, sigma0
from ase.calculators.qmmm import EIQMMM, LJInteractionsGeneral
from ase.calculators.turbomole import Turbomole
from ase.constraints import FixBondLengths
from ase.optimize import LBFGS
from ase.io import read
import ase.units as units
import pandas as pd
from ase.calculators.loggingcalc import LoggingCalculator # EP

f = open('qmmm-znwat-LJgen.txt','w')  # Added by EP

mol = read('wat-zn.xyz')

mm_idx = list(range(0,(len(mol) - 48)))
qm_idx = list(range((len(mol) - 48),len(mol)))

#parameters for mm atoms (water) from TIP3
sigmamm  = np.array([sigma0,0,0])
epsilonmm = np.array([epsilon0,0,0])

#from prm file - epsilon is the second column (+) sigma is the third column
#1) get a list of qm atom types..

qm_atom_types = ['CT1','CT2','CPH1','NR1','CPH1','CPH2','NR2','CT1','CT2','CPH1','NR1','CPH1','CPH2','NR2','CT1','CT2','CPH1','NR1','CPH1','CPH2','NR1','ZN','OH1','H','HA','HA','HR3','HR1','H','HA','HA','HR3','HR1','H','HA','HA','HR3','HR1','H','H','H','H','H','H','H','H','H','H']

#2) get values for qm_atom_types
df = pd.read_csv('charmm22_param', sep=' *', error_bad_lines=False, header=None, usecols=[0,2,3], engine='python')
df.columns = ['atom_type','min_epsilon','sigma']

epsilonqm = []
for atomtype in qm_atom_types:
    epsilonqm.append((-(df.loc[df['atom_type'] == atomtype, 'min_epsilon'].iloc[0]))* units.kcal / units.mol)
epsilonqm = np.array(epsilonqm)

sigmaqm = []
for atomtype in qm_atom_types:
    sigmaqm.append(((df.loc[df['atom_type'] == atomtype, 'sigma'].iloc[0])*2)*(2**(-1/6)))
sigmaqm = np.array(sigmaqm)

#insert values in arrays as interactions with Class LJInteractionsGeneral
interaction = LJInteractionsGeneral(sigmaqm, epsilonqm, epsilonmm, sigmamm, qm_molecule_size=48)

qm_par = {'title': 'zn cluster in wat qmmm'}

mol.set_pbc(True)

mol.calc =(EIQMMM(qm_idx, Turbomole(restart=True,**qm_par), TIP3P(), interaction, vacuum=4.0, output=f))

print(mol.get_potential_energy(), file = f)



mm_bonds = [(3 * i + j, 3 * i + (j + 1) % 3) for i in range(int((len(mol)-2) / 3)) for j in [0, 1, 2] if (3 * i + j) not in qm_idx]
mol.constraints = FixBondLengths(mm_bonds)

opt = LBFGS(mol, trajectory='zn_inwatbox_LJgen_1.traj', restart='zn_inwatbox_LJgen_1.pckl')
opt.run(fmax=0.05)


f.close()
