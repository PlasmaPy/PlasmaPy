"""
This module tests the functions in the maagnetics class
"""
import h5py
from matplotlib import pyplot as plt
from magnetics import Magnetics

def printname(name):
    """this function does sth"""
    print(name)

def load_hdf5(file,verbose=False):
    """this function also does sth"""
    f = h5py.File(file, 'r')
    if verbose:
        print('All Groups Contained')
        f.visit(printname)
    return f

directory = '/Users/aminaahmed/Documents/Project_Plasma/'
datafilename = '6202022_PracticeShots_9total_5or6areDecent.h5'
data = load_hdf5(directory + datafilename ,verbose=True)

test = Magnetics()
Bdot = data['pos19']['Bdot']['theta'][1,:]
tloop = 1.6129508E-6
times = data['time']['time_s'][0,:]

result = test.bdot_field(Bdot, tloop, times)

#plots
timeB_s = data['time']['timeB_s'][0,:]
comp = data['pos19']['B']['theta'][1,:]
plt.plot(timeB_s, comp)
plt.title('Magnetic Field of Shot')
plt.show()
plt.plot(timeB_s, result)
plt.title("Magnetic field calculated")
plt.show()
