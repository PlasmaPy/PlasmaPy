__all__ = [
    "AbstractClassicalTransportCoefficients",
]



import numpy as np
import astropy.units as u
from abc import ABC

from plasmapy.formulary.parameters import gyrofrequency
from plasmapy import particles
from plasmapy.particles import Particle



class AbstractClassicalTransportCoefficients(ABC):
    r"""
    An abstract class representing classical transport coefficients. 
    
    Subclasses representing different classical transport models re-implement
    the transport cofficient methods of this abstract class.
    """
    
    #@particles.particle_input
    def __init__(self, chi, Z,
                 B : u.T = None,
                 n: u.m**-3 = None,
                 T: u.K = None,
                 particle: Particle = '',
                 ):
           
        self.chi = chi
        self.Z = Z
        
        self.B = B
        self.n = n
        self.T = T
        self.particle = particle

        try:
            self.ei_collision_period =  chi / gyrofrequency(B, 'e-')
        except ValueError:
            self.ei_collision_period = None
        
        # Normalizations are defined such that multiplying the dimensionless
        # quantity by the normalization constant will return the dimensional
        # quantity
        try:
            self.alpha_normalization = self.particle.mass * self.n / self.ei_collision_period
        except (TypeError, AttributeError):
            self.alpha_normalization = None
            
        self.beta_normalization = 1.0 # Beta is naturally dimensionless
        
        try:
            self.kappa_normalization = self.n * self.T * self.ei_collision_period / self.particle.mass
        except (TypeError, AttributeError): 
            self.kappa_normalization = None

        
        
    # **********************************************************************
    # Resistivity (alpha)
    # **********************************************************************
    def norm_alpha_para(self):
        raise NotImplementedError
        
    def alpha_para(self):
        if self.alpha_normalization is None:
            raise ValueError("Keywords n, B, and particle must be provided to "
                             "calculate alpha_para")
            
        return self.norm_alpha_para(self.chi, self.Z) * self.alpha_normalization

    def norm_alpha_perp(self):
        raise NotImplementedError
        
    def alpha_perp(self):
        if self.alpha_normalization is None:
            raise ValueError("Keywords n, B, and particle must be provided to "
                             "calculate alpha_perp")
            
        return self.norm_alpha_perp(self.chi, self.Z) * self.alpha_normalization

    def norm_alpha_cross(self):
        raise NotImplementedError
        
    def alpha_cross(self):
        if self.alpha_normalization is None:
                    raise ValueError("Keywords n, B, and particle must be provided to "
                                     "calculate alpha_cross")
        return self.norm_alpha_cross(self.chi, self.Z) * self.alpha_normalization
    
    
    # **********************************************************************
    # Thermoelectric Coefficient (beta)
    # **********************************************************************
    
    
    # **********************************************************************
    # Thermal Conductivity (kappa)
    # **********************************************************************
    
    