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


validMathTextFonts = ['dejavusans',
                      'dejavuserif',
                      'cm',
                      'stix',
                      'stixsans',
                      'custom']

# setting latex rendered fonts to be same as regular fonts
try:
    matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'
    matplotlib.rcParams['mathtext.rm'] = 'DejaVu Serif'
    matplotlib.rcParams['mathtext.it'] = 'DejaVu Serif:italic'
    matplotlib.rcParams['mathtext.bf'] = 'DejaVu Serif:bold'
except:
    print("Couldn't load dejavuserif fonts for plot defaults."
          "Falling back to stix fonts.")
    try:
        matplotlib.rcParams['mathtext.fontset'] = 'stix'
        matplotlib.rcParams['font.family'] = 'STIXGeneral'
    except:
        print("Couldn't load stix fonts for plot defaults.")


# comparing fonts
plt.title(r"cm-3 $\rm cm^{-3}$")
plt.show()


