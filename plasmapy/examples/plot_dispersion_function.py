"""
The plasma dispersion function
==============================

Let's import some basics (and `PlasmaPy`!)
"""


import numpy as np
import matplotlib.pyplot as plt
import plasmapy


#######################################################################
help(plasmapy.mathematics.plasma_dispersion_func)


#######################################################################
# We'll now make some sample data to visualize the dispersion function:

x = np.linspace(-1, 1, 1000)
X, Y = np.meshgrid(x, x)
Z = X + 1j * Y
print(Z.shape)

#######################################################################
# Before we start plotting, let's make a visualization function first:


def plot_complex(X, Y, Z, N=50):
    fig, (real_axis, imag_axis) = plt.subplots(1, 2)
    real_axis.contourf(X, Y, Z.real, N)
    imag_axis.contourf(X, Y, Z.imag, N)
    real_axis.set_title("Real values")
    imag_axis.set_title("Imaginary values")
    for ax in [real_axis, imag_axis]:
        ax.set_xlabel("Real values")
        ax.set_ylabel("Imaginary values")
    fig.tight_layout()


plot_complex(X, Y, Z)

#######################################################################
# We can now apply our visualization function to our simple

F = plasmapy.mathematics.plasma_dispersion_func(Z)
plot_complex(X, Y, F)


#######################################################################
# So this is going to be a hack and I'm not 100% sure the dispersion function
# is quite what I think it is, but let's find the area where the dispersion
# function has a lesser than zero real part because I think it may be important
# (brb reading Fried and Conte):

plot_complex(X, Y, F.real < 0)


#######################################################################
# We can also visualize the derivative:

F = plasmapy.mathematics.plasma_dispersion_func_deriv(Z)
plot_complex(X, Y, F)

#######################################################################
# Plotting the same function on a larger area:

x = np.linspace(-2, 2, 2000)
X, Y = np.meshgrid(x, x)
Z = X + 1j * Y
print(Z.shape)

#######################################################################

F = plasmapy.mathematics.plasma_dispersion_func(Z)
plot_complex(X, Y, F, 100)

#######################################################################
# Now we examine the derivative of the dispersion function as a function
# of the phase velocity of an electromagnetic wave propagating through
# the plasma. This is recreating figure 5.1 in:
# J. Sheffield, D. Froula, S. H. Glenzer, and N. C. Luhmann Jr,
# Plasma scattering of electromagnetic radiation: theory and measurement
# techniques. Chapter 5 Pg 106 (Academic press, 2010).

xs = np.linspace(0, 4, 100)
ws = (-1 / 2) * plasmapy.mathematics.plasma_dispersion_func_deriv(xs)
wRe = np.real(ws)
wIm = np.imag(ws)

plt.plot(xs, wRe, label="Re")
plt.plot(xs, wIm, label="Im")
plt.axis([0, 4, -0.3, 1])
plt.legend(loc='upper right',
           frameon=False,
           labelspacing=0.001,
           fontsize=14,
           borderaxespad=0.1)
plt.show()