#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 14:16:48 2018

Helper functions for making plots

@author: Pawel M. Kozlowski
"""

import matplotlib.pyplot as plt

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