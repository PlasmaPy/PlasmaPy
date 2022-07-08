# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 17:01:56 2020

@author: Josh0
"""

#load_hdf5.py
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:44:36 2019
@author: dschaffner
"""

import h5py

def printname(name):
    print(name)

def load_hdf5(file,verbose=False):
    f = h5py.File(file, 'r')
    if verbose:
        print('All Groups Contained')
        f.visit(printname)
    return f