#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 14:01:17 2020

@author: elron
"""

''' script to get all distances and angles of the qm system'''


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


def distances(starttitle, endtitle, qm_idx):
    
    '''
    starttitle: string
        title of the initial structure e.g cluster in waterbox at t = 0
    
    endtitle: string
        title of the trajectory movement file of the cluster in waterbox
        
    qm_idx: integer
        number of qm atoms with negative algebraic sign
    '''
    

    mol0 = read(starttitle)[qm_idx:]
    mol1 = Trajectory(endtitle)[-1][qm_idx:]
    symbollist = list(mol0.symbols)
    # creates columnnames with atom+index e.g C1
    newlist = []
    counter = 0
    for i in symbollist:
        counter = (counter+1)
        a = str(counter)
        newlist.append(i+a)
    
    ### Dataframe for initial structure

    symbols = pd.DataFrame(list(mol0.symbols), columns = ['Atomtype'], index = [i for i in newlist])
    start_dist = pd.DataFrame(mol0.get_all_distances(), index = [i for i in newlist], columns = newlist)#.join(symbols)
    
     ### Dataframe for final structure

    symbols = pd.DataFrame(list(mol1.symbols), columns = ['Atomtype'], index = [i for i in newlist])
    end_dist = pd.DataFrame(mol1.get_all_distances(), index = [i for i in newlist], columns = newlist)#.join(symbols)
    

    # Difference between final and initial structure
    diff_dist = start_dist-end_dist
    diff_dist.to_csv('diff_dist.csv')
    start_dist.to_csv('start_dist.csv')
    end_dist.to_csv('end_dist.csv')
    
    #return start_dist, end_dist, diff_dist

    # Displacment from startstruktur in percent
    # formula x2-x1 = delta
    # x1/delta = Displacment bezogen auf Startstruktur
    
    percent_dist = (diff_dist/start_dist)*100
    percent_dist.to_csv('percent_dist.csv')
    
    
    
   
distances('znclus.xyz','zn_inwatbox_LJgen_1.traj',-48)

diff = pd.read_csv('test.csv').set_index("Unnamed: 0")
start = pd.read_csv('start.csv').set_index("Unnamed: 0")
end = pd.read_csv('end.csv').set_index("Unnamed: 0")
percent = pd.read_csv('percent.csv').set_index("Unnamed: 0")
#dfiff = df.set_index("Unnamed: 0")

