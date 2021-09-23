import os
import numpy as np

from plasmapy.data import data

def test_get_file():
    filename = 'NIST_PSTAR_aluminum.txt'
    
    # Download data (or check that it already exists)
    path = data.get_file(filename)
    
    # verify contents match expectations
    arr = np.loadtxt(path, skiprows=7)
    assert np.allclose(arr[0,:], np.array([1e-3, 1.043e2]))
    
    # delete the file
    os.remove(path)
    
    # Re-download
    path = data.get_file(filename)
