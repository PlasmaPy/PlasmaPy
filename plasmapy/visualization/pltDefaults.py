#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 02:37:12 2017

Default plotting parameters

@author: Pawel M. Kozlowski
"""

# importing Python modules
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from itertools import cycle

# setting default plot properties
plotFont = {'family' : 'serif',
            'weight' : 'normal',
            'size'   : 14}
plt.rc('font', **plotFont)
plt.rc('lines', linewidth=2)

# setting latex rendered fonts to be same as regular fonts
matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'
matplotlib.rcParams['mathtext.rm'] = 'DejaVu Serif'
matplotlib.rcParams['mathtext.it'] = 'DejaVu Serif:italic'
matplotlib.rcParams['mathtext.bf'] = 'DejaVu Serif:bold'

# for testing whether matplotlib and python fonts match
#plt.title(r'cm$\rm cm^{-3}$')


#%% convenience functions for plotting

# shaded error bars for line plot
def plot_line_shaded(xData, yData, yErrs, label="", **kwargs):
    """
    Generate a line plot with shaded region representing y-error bars.
    Can be run multiple times before plt.show(), to plot multiple data
    sets on the same axes.
    x axis data points
    y axis data points
    y axis errors
    """
    # check that arrays are 1D
    
    # check that arrays are equal length
    if not len(xData) == len(yData) == len(yErrs):
        raise ValueError("Arrays must of equal length!")
    
    # line plot
    ax = kwargs.pop('ax', plt.gca())
    base_line, = ax.plot(xData, yData, label=label, **kwargs)
    # shaded error region
    plt.fill_between(xData,
                     yData - yErrs,
                     yData + yErrs,
                     alpha=0.5,
                     edgecolor=base_line.get_color(),
                     facecolor=base_line.get_color())
    
    return

# error bars for scatter plot
def plot_scatter_bars(xData, yData, yErrs, label="", **kwargs):
    """
    Generate a scatter plot with y-error bars.
    Can be run multiple times before plt.show(), to plot multiple data
    sets on the same axes.
    x axis data points
    y axis data points
    y axis errors
    """
    # check that arrays are 1D
    
    # check that arrays are equal length
    if not len(xData) == len(yData) == len(yErrs):
        raise ValueError("Arrays must of equal length!")
        
    # line plot
    ax = kwargs.pop('ax', plt.gca())
    ax.errorbar(xData,
                yData,
                yerr=yErrs,
                label=label,
                fmt='o',
                fillstyle='none',
                capsize=4,
                elinewidth=1,
                **kwargs)
    return

#%% testing custom plotting tools

#xData = np.arange(10)
#yData = xData ** 2
#yData2 = xData ** 3
#yErrs = np.ones_like(xData) * 30
#plot_line_shaded(xData, yData, yErrs, label="#1", color='r')
#plot_line_shaded(xData, yData2, yErrs, label="#2")
#plt.ylabel("y")
#plt.xlabel("x")
#plt.legend(loc='upper left', frameon=False, labelspacing=0.001,
#           fontsize=14, borderaxespad=0.4)
#plt.show()
#
#
#xData = np.arange(10)
#yData = xData ** 2
#yData2 = xData ** 3
#yErrs = np.ones_like(xData) * 30
#plot_scatter_bars(xData, yData, yErrs, label="#1", color='r')
#plot_scatter_bars(xData, yData2, yErrs, label="#2")
#plt.ylabel("y")
#plt.xlabel("x")
#plt.legend(loc='upper left', frameon=False, labelspacing=0.001,
#           fontsize=14, borderaxespad=0.4)
#plt.show()

#%% style cycling
# printing default colors
default_colors = matplotlib.colors.cnames.keys()
#print(f"Colors: {default_colors}")
# printing default linestyles
default_lines = matplotlib.lines.lineStyles.keys()
#print(f"Lines: {default_lines}")
# printing default marker styles
default_markers = matplotlib.markers.MarkerStyle.markers.keys()
#print(f"Markers: {default_markers}")

# custom list of linestyles (excluding blank line styles)
lines = ["-","--","-.",":"]
linecycler = cycle(lines)
#plt.figure()
#for i in range(10):
#    x = range(i,i+10)
#    plt.plot(range(10),x,next(linecycler))
#plt.show()
    