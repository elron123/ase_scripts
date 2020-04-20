#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 17:18:27 2020

@author: elron
"""

import numpy as np
from ase.calculators.tip3p import TIP3P, epsilon0, sigma0
from ase.calculators.qmmm import EIQMMM, LJInteractionsGeneral
from ase.calculators.turbomole import Turbomole
from ase.constraints import FixBondLengths
from ase.optimize import LBFGS
from ase.io import read
import ase.units as units
import pandas as pd
from ase.io import Trajectory, read, write
from ase.visualize import view
import re
import os
import PyPDF2
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileMerger, PdfFileReader

#### Code snippet for graph making

import matplotlib.pyplot as plt
import seaborn as sns

bp86 = pd.read_csv('znclus_bp86.csv') 
pbe = pd.read_csv('znclus_pbe.csv')

bp86 = bp86[bp86['intensity_%']!=0]
pbe = pbe[pbe['intensity_%']!=0]

# def ir_graph(title, df):

    
    
    # plt.figure(figsize = [15,10])
    # plt.plot(pbe['frequency'],pbe['intensity_%'], color = 'k')
    # plt.title('IR spectra {}'.format(title), fontsize = 18)
    # plt.xlabel('Wavenumbers [cm-1]', fontsize = 12)
    # plt.ylabel('Intensity [%]', fontsize = 12)
    # ax.set_xlim(0,3670)
    # ax.set_xticklabels([str(i) for i in ax.get_xticks()],fontsize = 12)
    # ax.set_yticklabels([str(i) for i in ax.get_xticks()],fontsize = 12)
    # ax.invert_yaxis()
    # ax.invert_xaxis()
    # fig = ax.get_figure()
    # fig.savefig('{}.png'.format(title))
    # plt.show()





fig, ax = plt.subplots(figsize=(15,10))
color = 'tab:red'
ax.plot(bp86['frequency'],bp86['intensity_%'], color=color, label ='[Zn]-BP86')
plt.title('Overlapped IR spectra of the Zncluster' ,fontsize = 18)
plt.xlabel('Wavenumbers [cm-1]', fontsize = 12)
plt.ylabel('Intensity [%]', fontsize = 12)
ax.legend(loc = 'lower left')
ax.invert_yaxis()
#ax.invert_xaxis()

#ax = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax.plot(pbe['frequency'],pbe['intensity_%'], color=color,label ='[Zn]-PBE')
#ax.invert_yaxis()
#ax.invert_xaxis()
fig.tight_layout()
ax.legend(loc = 'lower left')
#g = fig.get_figure()
fig.savefig('test.png')
plt.show()






