import math
import numpy as np
import scipy.special
import sympy.functions.special.bessel

from sympy import Derivative  
  
class ForceFreeFluxRope:
  
    def __init__(self, B0, a, lamb):
        self.B0 = B0
        self.a = a
        self.lamb = lamb
        self.r = r
        self.z = z



    def B_radial_cyln():

            return (0)

    def B_phi_cyln(self, B, r):

            return (B*scipy.special.j1(self.lamb * r))

    def B_z_cyln(self, B, r):

            return(B*scipy.special.j1(self.lamb * r))