"""
This module tests the functions in the maagnetics class
It is a scratch file
will be removed 
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
timeB_s = data['time']['timeB_s'][0,:]
time_s = data['time']['time_s'][0,:]

"""
#bdot_field function 
Bdot = data['pos19']['Bdot']['theta'][1,:]
tloop = 1.6129508E-6
times = data['time']['time_s'][0,:]

result = test.bdot_field(Bdot, tloop, times, "Tesla")

#plots
timeB_s = data['time']['timeB_s'][0,:]
comp = data['pos19']['B']['theta'][1,:]

plt.plot(timeB_s, comp)
plt.title('Magnetic Field of Shot')
plt.show()

plt.plot(timeB_s, result)
plt.title("Magnetic field calculated")
plt.show()
"""

#finding screwed up shots
#pos19
b19_r = data['pos19']['B']['r']

"""
count = 0
for arr in b19_r:
    if count == 0: 
        count = count + 1
        continue
    plt.plot(timeB_s, arr)
    plt.title("Array number: " + str(count))
    plt.show()  
    count = count + 1
"""

test2_array = b19_r[1,:]
subbed_array = test.band_pass_filter(test2_array, 5e4, 125e6, 3)

plt.plot(timeB_s, test2_array)
plt.title("Array before filter" )
plt.show() 

plt.plot(timeB_s, subbed_array)
plt.title("Array after filter" )
plt.show()   

comp_data = data['pos19']['Bdot']['r'][1,:]
plt.plot(comp_data)
plt.title("coresponding bdot data" )
plt.show()  