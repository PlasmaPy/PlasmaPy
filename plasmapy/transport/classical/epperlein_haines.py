"""
Implementation of the magnetized transport coefficients presented by
Epperlein and Haines in their 1986 paper "Plasma transport coefficients in a magnetic field by direct numerical solution of
the Fokker–Planck equation" (Phys. Fluids), doi: 10.1063/1.865901
"""

__all__ = [
    "EpperleinHainesCoefficents",
    ]


# TODO: allow options for returning transport coefficients in physical units
# this will require more information about the plasma
# look at sec. IIB in Ridgers2008 (pg 3) for relations: these likely hold for EH as well
# EH equivalent is right at the end of section II (pg. 5). Looks the same.

from abc import ABC
import astropy.units as u
import numpy as np
import warnings

from plasmapy.formulary.parameters import gyrofrequency
from plasmapy import particles
from plasmapy.particles import Particle

from plasmapy.transport.classical.base import AbstractClassicalTransportCoefficients

import matplotlib.pyplot as plt

Ztable = np.array([1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 20, 30, 60, np.inf])

# TABLE II
α0 = np.array([0.5061, 0.4295, 0.3950, 0.3750, 0.3618, 0.3524, 0.3454, 0.3399, 0.3319, 0.3263, 0.3221, 0.3144, 0.3081, 0.3015, 0.2945])
α0p = np.array([1.37, 1.58, 1.68, 1.74, 1.78, 1.8, 1.82, 1.84, 1.87, 1.88, 1.9, 3.53, 5.49, 7.61, 9.17])
α1p = np.array([3.03, 3.21, 3.17, 3.15, 3.14, 3.13, 3.12, 3.11, 3.1, 3.1, 3.09, 3.52, 3.97, 4.41, 4.73])
a0p = np.array([2.77, 2.78, 2.78, 2.78, 2.78, 2.79, 2.79, 2.79, 2.79, 2.8, 2.8, 5.14, 7.94, 10.9, 13])
a1p = np.array([6.72, 6.7, 6.47, 6.37, 6.33, 6.29, 6.26, 6.23, 6.21, 6.2, 6.19, 7.97, 10, 12.2, 13.8])
α0x = np.array([0.1989, 0.3285, 0.41457, 0.4794, 0.5285, 0.5676, 0.5996, 0.6264, 0.6686, 0.7004, 0.7254, 0.7757, 0.8209, 0.8727, 0.9328])
α0pp = np.array([2.66e2, 4.91e2, 6.3e2, 7.06e2, 7.23e2, 7.57e2, 7.94e2, 8.17e2, 8.89e2, 9.62e2, 9.88e2, 1.06e3, 1.17e3, 1.3e3, 1.33e3])
α1pp = np.array([2.53] * 15)
a0pp = np.array([3280, 3720, 3790, 3670, 3370, 3280, 3250, 3200, 3270, 3390, 3360, 3360, 3520, 3720, 3540])
a1pp = np.array([3460, 6690, 8990, 10200, 10600, 11100, 11700, 12200, 13300, 14700, 15100, 16200, 18500, 20800, 21900])
a2pp = np.array([366, 554, 675, 735, 744, 769, 793, 825, 863, 925, 940, 985, 1080, 1180, 1180])

# TABLE III
#β0(Z=20) is clearly incorrect; replaced with fitted value
β0 = np.array([0.702, 0.9054, 1.018, 1.092, 1.146, 1.186, 1.218, 1.244, 1.283, 1.312, 1.334, 1.3801, 1.414, 1.455, 1.5])
β0p = np.array([1.05e3, 1.38e3, 1.55e3, 1.64e3, 1.71e3, 1.74e3, 1.73e3, 1.79e3, 1.92e3, 1.89e3, 1.92e3, 2.01e3, 2.09e3, 2.16e3, 2.20e3])
β1p = np.array([6.33] * 15)
b0p = np.array([3710, 3800, 3800, 3740, 3720, 3640, 3530, 3570, 3740, 3580, 3580, 3620, 3680, 3690, 3650])
b1p = np.array([4110, 7050, 8750, 9730, 10700, 11100, 11300, 11700, 12900, 13000, 13300, 14300, 15200, 16000, 16800])
b2p = np.array([515, 642, 690, 707, 731, 735, 729, 734, 772, 762, 765, 789, 808, 822, 809])
β0x = np.array([0.8831, 1.812, 2.589, 3.236, 3.78, 4.243, 4.64, 4.986, 5.556, 6.007, 6.372, 7.144, 7.874, 8.756, 9.844])
β0pp = np.array([2.54, 4.40, 3.77, 3.43, 3.20, 3.05, 2.92, 2.82, 2.70, 2.62, 2.55, 2.43, 2.31, 2.26, 2.15])
β1pp = np.array([1.5] * 15)
b0pp = np.array([2.87, 2.43, 1.46, 1.06, .848, .718, .629, .565, .486, .436, .401, .341, .294, .258, .219])
b1pp = np.array([3.27, 5.18, 4.34, 3.92, 3.66, 3.48, 3.33, 3.21, 3.09, 3.00, 2.93, 2.81, 2.68, 2.64, 2.53])
b2pp = np.array([7.09, 9.34, 8.65, 8.27, 8.02, 7.83, 7.68, 7.55, 7.41, 7.30, 7.22, 7.07, 6.91, 6.84, 6.72])


# TABLE IV
γ0 = np.array([3.203, 4.931, 6.115, 6.995, 7.68, 8.231, 8.685, 9.067, 9.673, 10.13, 10.5, 11.23, 11.9, 12.67, 13.58])
γ0p = np.array([6.18, 9.30, 10.2, 9.14, 8.60, 8.57, 8.84, 7.93, 7.44, 7.32, 7.08, 6.79, 6.74, 6.36, 6.21])
γ1p = np.array([4.66, 3.96, 3.72, 3.6, 3.53, 3.49, 3.49, 3.43, 3.39, 3.37, 3.35, 3.32, 3.3, 3.27, 3.25])
c0p=np.array([1.93, 1.89, 1.66, 1.31, 1.12, 1.04, 1.02, 0.875, 0.77, 0.722, 0.674, 0.605, 0.566, 0.502, 0.457])
c1p=np.array([2.31, 3.78, 4.76, 4.63, 4.62, 4.83, 5.19, 4.74, 4.63, 4.7, 4.64, 4.65, 4.81, 4.71, 4.81])
c2p= np.array([5.3, 7.78, 8.88, 8.8, 8.8, 8.96, 9.24, 8.84, 8.71, 8.73, 8.65, 8.6, 8.66, 8.52, 8.53])
γ0x = np.array([6.071, 15.75, 25.65, 34.95, 43.45, 51.12, 58.05, 64.29, 75.04, 83.93, 91.38, 107.8, 124.3, 145.2, 172.7])
γ0pp = np.array([4.01, 2.46, 1.13, 0.628, 0.418, 0.319, 0.268, 0.238, 0.225, 0.212, 0.202, 0.200, 0.194, 0.189, 0.186])
γ1pp = np.array([2.5] * 15)
c0pp = np.array([.661, .156, .0442, .018, .00963, .00625, .00461, .00371, .003, .00252, .00221, .00185, .00156, .0013, .00108])
c1pp = np.array([.931, .398, .175, .101, .0702, .0551, .0465, .041, .0354, .0317, .0291, .0256, .0228, .0202, .018])
c2pp = np.array([2.5, 1.71, 1.05, .775, .646, .578, .539, .515, .497, .482, .471, .461, .45, .44, .43])


def _find_nearest(Z):
    """
    Finds the nearest Z-value in the EH tables for a given Z. Prints a 
    warning if the Z found is not equal to the Z requested.

    Parameters
    ----------
    Z : float
        An integer charge

    Returns
    -------
    i : int
        The index of the closest Z in the EH tables

    """
    i = np.argmin(np.abs(Ztable - Z))
    if Ztable[i] != Z:
        warnings.warn(f"Value Z = {Z} is not in the Epperlein-Haines table. "
                      f"Using the nearest value, Z = {Ztable[i]}. "
                      f"The values in the table are {Ztable}.",
                      RuntimeWarning)
    return i


class EpperleinHainesCoefficents(AbstractClassicalTransportCoefficients):
    
    def norm_alpha_para(self):
        """
        Calculates the normalized alpha_para coefficient in terms of the
        dimensionless Hall parameter and the ionization fraction.
    
        Parameters
        ----------
        chi : float 
            The dimensionless hall parameter (ratio of the electron gyrofrequency
            and the electron-ion collision frequency).
            
        Z : float
            Ionization fraction. The value will be coerced to the nearest value
            from the Epperlein-Haines tables: 
            Z = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 20, 30, 60, ∞]
    
        Returns
        -------
        norm_alpha_para, u.Quantity
            The resistivity parallel to (or in the absence of) a magnetic field
            
        Notes
        -----
        For details, see "Plasma transport coefficients in a magnetic field by 
        direct numerical solution of the Fokker–Planck equation" by
        E. M. Epperlein and M. G. Haines, DOI: `10.1063/1.865901`
    
        .. _`10.1063/1.865901`: https://aip.scitation.org/doi/10.1063/1.865901
    
        """
        i = _find_nearest(Z)
        return α0[i]*np.ones(self.chi.size)


    def norm_alpha_perp(self):
        i = _find_nearest(self.Z)
        return 1 - (α1p[i]*self.chi + α0p[i])/(self.chi**2 + a1p[i]*self.chi + a0p[i])
    
    def norm_alpha_cross(self):
        i = _find_nearest(self.Z)
        return self.chi*(α1pp[i]*self.chi + α0pp[i])/ \
                    (self.chi**3 + a2pp[i]*self.chi**2 + a1pp[i]*self.chi + a0pp[i])**(8/9)
    
    def norm_beta_para(self):
        i = _find_nearest(self.Z)
        return β0[i]*np.ones(self.chi.size)
    
    def norm_beta_perp(self):
        i = _find_nearest(self.Z)
        return (β1p[i]*self.chi + β0p[i]) / (self.chi**3 + b2p[i]*self.chi**2 + b1p[i]*self.chi + b0p[i] )**(8/9)
    
    # TODO: Note that this function doesn't match the EH paper Fig. 1 in the
    # chi -> inf side. The coefficients and polynomial are right...
    # this might be a mistake in the EH tables?
    def norm_beta_cross(self):
        i = _find_nearest(self.Z)
        return self.chi*(β1pp[i]*self.chi + β0pp[i]) /  \
                    (self.chi**3 * b2pp[i]*self.chi**2 + b1pp[i]*self.chi + b0pp[i])
                    
    def norm_kappa_para(self):
        i = _find_nearest(self.Z)
        return γ0[i]*np.ones(self.chi.size)
                    
    def norm_kappa_perp(self):
        i = _find_nearest(self.Z)
        return (γ1p[i]*self.chi + γ0p[i])/ \
                    (self.chi**3 + c2p[i]*self.chi**2 + c1p[i]*self.chi + c0p[i])
    
    def norm_kappa_cross(self):
        i = _find_nearest(self.Z)
        return self.chi*(γ1pp[i]*self.chi + γ0pp[i])/  \
                 (self.chi**3 + c2pp[i]*self.chi**2 + c1pp[i]*self.chi + c0pp[i])






if __name__ == '__main__':
    chi = np.linspace(-2, 2, num=100)
    chi = 10**chi

    
    # Instantiate the object
    coef1 = EpperleinHainesCoefficents(chi, 1)
    coef2 = EpperleinHainesCoefficents(chi, np.inf)
    
   

    fig, axarr = plt.subplots(nrows=3, ncols=2, figsize=(10,10), sharex=True)
    
    for i in range (3):
        for j in range(2):
            ax = axarr[i][j]
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(1e-2, 1e2)
        
    axarr[-1][0].set_xlabel("$\\chi$")
    axarr[-1][1].set_xlabel("$\\chi$")
    
    ax = axarr[0][0]
    ax.set_title('$\\alpha_\\perp$')
    ax.set_ylim(0.2, 1)
    ax.plot(chi, coef1.norm_alpha_perp(), label="Z=1")
    ax.plot(chi, coef2.norm_alpha_perp(),  label="Z=Infinite")
    ax.legend()
    
    ax = axarr[1][0]
    ax.set_title('$\\alpha_\\wedge$')
    ax.set_ylim(0.01, 1)
    ax.plot(chi, coef1.norm_alpha_cross())
    ax.plot(chi, coef2.norm_alpha_cross())
    
    ax = axarr[2][0]
    ax.set_title('$\\beta_\\perp$')
    ax.set_ylim(0.001, 10)
    ax.plot(chi, coef1.norm_beta_perp())
    ax.plot(chi, coef2.norm_beta_perp())
    
    
    ax = axarr[0][1]
    ax.set_title('$\\beta_\\wedge$')
    ax.set_ylim(0.01, 1)
    ax.plot(chi, coef1.norm_beta_cross())
    ax.plot(chi, coef2.norm_beta_cross())
    
    ax = axarr[1][1]
    ax.set_title('$\\kappa_\\perp$')
    ax.set_ylim(0.001, 100)
    ax.plot(chi, coef1.norm_kappa_perp())
    ax.plot(chi, coef2.norm_kappa_perp())
    
    ax = axarr[2][1]
    ax.set_title('$\\kappa_\\wedge$')
    ax.set_ylim(0.01, 10)
    ax.plot(chi, coef1.norm_kappa_cross())
    ax.plot(chi, coef2.norm_kappa_cross())
    
