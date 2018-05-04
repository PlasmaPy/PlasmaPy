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
