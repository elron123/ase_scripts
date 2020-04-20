#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 10:13:41 2020

@author: elron
"""

########## Function combine water and cluster

# Needed Libaries
from ase.io import read, Trajectory, write
from ase.visualize import view

def combine_water_cluster(trajfile, clusterfile,outputfile, dist_qm_mm = 2):
    
    '''
        trajfile: string
            Trajectory of the atomic movement of the waterbox in .traj format
        
        clusterfile: string
            Clustercoordinates in .xyz format
            
        dist_qm_mm: integer or float 
            distance to qm atoms, 
            per default set to 2
            
        outputfile: string  
            outputfile with cluster in waterbox in xyz.format
    '''

    ### Waterbox
    #import traj
    traj = Trajectory(trajfile)
    #set watbox to last point of traj
    watbox = traj[-1]
    #get cell of watbox
    Wcell = watbox.get_cell()
    
    ### Cluster
    #import cluster coordinates
    qm = read(clusterfile)
    #set cluster cell to watbox value
    qm.set_cell(Wcell)
    #center cluster in the cell
    qm.center(vacuum=Wcell[0][0]/3)
    
    ###merge watbox and cluster
    all = watbox.copy()
    all.extend(qm)
    all.write('temp.xyz')
    
    #define index of atoms in mm part and index of atoms in qm part
    mm_idx = list(range(0,(len(all) - len(qm))))
    qm_idx = list(range((len(all) - len(qm)),len(all)))
    
    #get distances between atoms in mm and atoms in qm
    #get index of atoms within 2Ang from any qm atom
    dist = []
    conn_list = []
    
    for i in mm_idx:
        for j in qm_idx:
            dist.append(all.get_distance(i, j))
            
    for i, distance in enumerate(dist):
        i = int(i/len(qm_idx))  # for every i-nth 
        if distance <=dist_qm_mm:
            conn_list.append(i)
    
    conn_set = set(conn_list)
    conn_list = list(conn_set)
        
    #get index of oxygens within set distance in Ang from any qm atom
    oindex = []
    
    for i in conn_list:
        if all[i].symbol == 'O':
            oindex.append(i)
            
    #get index of water molecules within set distance in Ang from any qm qtom
    watslist = []
    for oxy in oindex:
        watslist.append(oxy)
        watslist.append(oxy +1)
        watslist.append(oxy +2)
    
    watsset = set(watslist)
    watslist = list(watsset)
    
    #delete water molecules within 2Ang from any qm atom from all
    del all.constraints
    del all[watslist]

    #write all to a file
    all.write(outputfile+".xyz")
    
    
    
    
    
    