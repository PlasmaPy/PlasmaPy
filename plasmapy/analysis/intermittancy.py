"""
Created on Tue Aug  2 22:18:16 2022

@author: dschaffner
"""
import numpy as np

def generate_deltas_list(timeseries,dt,timestep):
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries)
    initial_index = 0
    final_index = int(indexstep)
    initial = timeseries[initial_index]
    final = timeseries[final_index]
    numloops = 0
    while initial_index < (total-(indexstep+1)):
        numloops+=1
        delta = final-initial
        deltas.append(delta)
        initial_index+=1
        final_index+=1
        initial = timeseries[initial_index]
        final = timeseries[final_index]
    return deltas