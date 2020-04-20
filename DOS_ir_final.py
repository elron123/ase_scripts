#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
@author: elron
"""

import numpy as np
import pandas as pd
import re
import os
from math import pi, sqrt, log

def make_ir_file(csvname):
    
    ### needed fixed names
    bad_words = ['mol']
    textnames = ['ir_extract.txt','ir_2.txt','ir_output.txt']
    TAG = 'WARNING: values of IR intensities are questionable for saddle point structures'
    
    ### Create first file and get all filterlist =['mode','frequency','symmetry','IR','intensity'] lines of aoforce output
    with open('ASE.TM.aoforce.out','r') as f, open('ir_extract.txt','w') as out:
        lines = f.readlines()
        for line in lines:
            if re.search(r'mode',line):
                print(' '.join(line.split()), file = out)   # join split formulation to remove empty lines and make just one blankspace between lines
            if re.search(r'frequency',line):
                print(' '.join(line.split()), file = out)
            if re.search(r'symmetry',line):
                print(' '.join(line.split()), file = out)
            if re.search(r'IR',line):
                print(' '.join(line.split()), file = out)  
            if re.search(r'intensity [( ]',line):
                line = line.replace('intensity (  %   )','intensity_%')
                print(' '.join(line.split()), file = out)
    
    ### Remove intensity (km/mol) lines with help of badword 'mol'
    with open('ir_extract.txt', 'r') as f, open('ir_2.txt','w') as out:
        for line in f:
            if not any(bad_word in line for bad_word in bad_words):
                out.write(line)

    ### remove unessarcy text of aoforce output file
    tag_found = False
    with open('ir_2.txt') as in_file, open('ir_output.txt', 'w') as out_file :
        for line in in_file:
            if not tag_found:
                if line.strip() == TAG:
                    tag_found = True
            else:
                out_file.write(line)
                
    ############### Create readable csv file from above made txt file
    
    ### make dataframe
    df = (pd.read_csv('ir_output.txt', sep = ' ', header = None, index_col = 0)
      .reset_index()
      .transpose()
      )

    df.columns = df.iloc[0]   # make first row column names
    df.drop(df.index[0], inplace = True) # drop first row
    filterlist =['mode','frequency','symmetry','IR','intensity_%'] 
    
    # append values of muliple same named columns, into unique named columns
    s = pd.Series(df.columns)
    df.columns = df.columns+s.groupby(s).cumcount().astype(str) # make columnnames unique
    df2 = (pd.wide_to_long(df.reset_index(),stubnames=filterlist,i='index',j='drop',suffix='\d+') #stubnames names of original columnsnames without attached number  
        .reset_index()
        .drop(columns = ['index','drop'], axis = 1)
        )
    
    # remove imaginary frequencies    
    df2[['frequency','intensity_%']] = df2[['frequency','intensity_%']].apply(pd.to_numeric, errors = 'coerce')
    df2.dropna(subset = ['frequency'], inplace = True)

    df2.to_csv(csvname+".csv")
    
    # delete all tempoary made files   
    for i in textnames:
        os.remove(i)



def fold(irfile,title,start=800.0, end=4000.0, npts=None, width=4.0,
             type='Gaussian', normalize=False):
        """Fold frequencies and intensities within the given range
        and folding method (Gaussian/Lorentzian).
        The energy unit is cm^-1.
        normalize=True ensures the integral over the peaks to give the
        intensity.
        """
        
        df = pd.read_csv(irfile)
        df = df[(df != 0).all(1)]
        intensities = df['intensity_%']
        with open('out.csv','w') as fd:        

            lctype = type.lower()
            assert lctype in ['gaussian', 'lorentzian']
            if not npts:
                npts = int((end - start) / width * 10 + 1)
            prefactor = 1
            if lctype == 'lorentzian':
                df['intensity_%'] = df['intensity_%'] * width * pi / 2.
                intensities = df['intensity_%'].copy()
                if normalize:
                    prefactor = 2. / width / pi
            else:
                sigma = width / 2. / sqrt(2. * log(2.))
                if normalize:
                    prefactor = 1. / sigma / sqrt(2 * pi)
        
            # Make array with spectrum data
            spectrum = np.empty(npts)
            energies = np.linspace(start, end, npts)
            for i, energy in enumerate(energies):
                energies[i] = energy
                if lctype == 'lorentzian':
                    spectrum[i] = (intensities * 0.5 * width / pi /
                                   ((df['frequency'] - energy)**2 +
                                    0.25 * width**2)).sum()
                else:
                    spectrum[i] = (intensities *
                                   np.exp(-(df['frequency'] - energy)**2 /
                                          2. / sigma**2)).sum()
            #return [energies, prefactor * spectrum]
    
            outdata = np.empty([len(energies), 2])
            outdata.T[0] = energies
            outdata.T[1] = spectrum
            #fd.write('# %s folded, width=%g cm^-1' % (type.title(), width))
            #fd.write('# [cm^-1] arbitrary')
            for row in outdata:
                fd.write('%.3f %9.5e\n' %
                         (row[0], row[1]))

            df1= pd.read_csv('out.csv',decimal = ",", sep=" ")
            df1.columns = ['frequency_cm^-1','dos']
            df1['frequency_cm^-1'] = df1['frequency_cm^-1'].astype(float)
            df1['frequency_cm^-1']=df1['frequency_cm^-1']
            df1['dos'] = df1['dos'].astype(float)
            df1.to_csv(title+"ir_spectra.csv")