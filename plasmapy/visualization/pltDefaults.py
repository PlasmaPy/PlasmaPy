"""Dynamically set defaults for matplotlib plots"""

# importing Python modules
import matplotlib 
import matplotlib.pyplot as plt

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