"""Tests for functions that calculate transport coefficients."""


import numpy as np
import pytest
from astropy import units as u
from plasmapy.atomic.atomic import ion_mass, charge_state

from ...utils.exceptions import RelativityWarning, RelativityError
from ...constants import c, m_p, m_e, e, mu0

from ..transport import (Coulomb_logarithm, classical_transport,
                         _nondim_tc_e_braginskii,
                         _nondim_tc_i_braginskii,
                         _nondim_tec_braginskii,
                         _nondim_resist_braginskii,
                         _nondim_visc_i_braginskii,
                         _nondim_visc_e_braginskii,
                         )

def count_decimal_places(digits):
    '''Return the number of decimal places of the input digit string'''
    integral, _, fractional = digits.partition(".")
    return len(fractional)


def test_Coulomb_logarithm():

    n_e = np.array([1e9, 1e9, 1e24])*u.cm**-3
    T = np.array([1e2, 1e7, 1e8])*u.K
    Lambda = np.array([5.97, 21.66, 6.69])
    particles = ('e', 'p')

    for i in range(3):
        assert np.isclose(Coulomb_logarithm(T[i], n_e[i], particles),
                          Lambda[i], atol=0.01)

    assert np.isclose(Coulomb_logarithm(1*u.eV, 5*u.m**-3, ('e', 'e')),
                      Coulomb_logarithm(11604.5220*u.K, 5*u.m**-3, ('e', 'e')))

    assert np.isclose(Coulomb_logarithm(1e2*u.K, 1e9*u.cm**-3, ('e', 'p')),
                      5.97, atol=0.01)

    assert np.isclose(Coulomb_logarithm(1e7*u.K, 1e9*u.cm**-3, ('e', 'p')),
                      21.6, atol=0.1)

    assert np.isclose(Coulomb_logarithm(1e8*u.K, 1e24*u.cm**-3, ('e', 'p')),
                      6.69, atol=0.01)

    assert np.allclose(Coulomb_logarithm(T, n_e, particles), Lambda, atol=0.01)

    assert np.isclose(Coulomb_logarithm(1e5*u.K, 5*u.m**-3, ('e', 'e'),
                                        V=1e4*u.m/u.s), 21.379082011)

    with pytest.warns(RelativityWarning):
        Coulomb_logarithm(1e5*u.K, 1*u.m**-3, ('e', 'p'), 0.9*c)

    with pytest.raises(RelativityError):
        Coulomb_logarithm(1e5*u.K, 1*u.m**-3, ('e', 'p'), 1.0*c)

    with pytest.raises(u.UnitConversionError):
        Coulomb_logarithm(1e5*u.g, 1*u.m**-3, ('e', 'p'), 29979245*u.m/u.s)

    with pytest.raises(ValueError):
        Coulomb_logarithm(1*u.K, 5*u.m**-3, ('e'))

    with pytest.raises(ValueError):
        Coulomb_logarithm(1*u.K, 5*u.m**-3, ('e', 'g'))

    with pytest.raises(ValueError):
        Coulomb_logarithm(1*u.K, 5*u.m**-3, ('e', 'D'))
        
        
# test class for classical_transport class:
class Test_classical_transport(object):

    def setup_method(self):
        """set up some initial values for tests"""
        u.set_enabled_equivalencies(u.temperature_energy())
        self.T_e = 1000 * u.eV
        self.n_e = 2e13 / u.cm ** 3
        self.ion_particle = 'D +1'
        self.m_i = ion_mass(self.ion_particle)
        self.Z = charge_state(self.ion_particle)
        self.T_i = self.T_e
        self.n_i = self.n_e / self.Z
        self.B = 0.01 * u.T
        self.coulomb_log_val_ei = 17
        self.coulomb_log_val_ii = 17
        self.hall_e = None
        self.hall_i = None
        self.V_ei = None
        self.V_ii = None
        self.mu = m_e / self.m_i
        self.theta = self.T_e / self.T_i
        self.model = 'Braginskii'
        self.field_orientation = 'all'
        self.ct = classical_transport(
            T_e = self.T_e, 
            n_e = self.n_e, 
            T_i = self.T_i, 
            n_i = self.n_i,
            ion_particle = self.ion_particle, 
            Z = self.Z,
            B = self.B, 
            model = self.model,
            field_orientation = self.field_orientation,
            coulomb_log_ei = self.coulomb_log_val_ei,
            coulomb_log_ii = self.coulomb_log_val_ii,
            V_ei = self.V_ei, 
            V_ii = self.V_ii,
            hall_e = self.hall_e, 
            hall_i = self.hall_i,
            mu = self.mu, 
            theta = self.theta,
        )

    def test_resistivity_units(self):
        """output should be a Quantity with units of Ohm m"""
        assert(self.ct.resistivity().unit == u.Ohm * u.m)
        
    def test_thermoelectric_conductivity_units(self):
        """output should be a Quantity with units of dimensionless"""
        assert(self.ct.thermoelectric_conductivity().unit == u.m / u.m)
        
    def test_ion_thermal_conductivity_units(self):
        """output should be Quantity with units of W / (m K)"""
        assert(self.ct.ion_thermal_conductivity().unit == u.W / u.m / u.K)
    
    def test_electron_thermal_conductivity_units(self):
        """output should be Quantity with units of W / (m K)"""
        assert(self.ct.electron_thermal_conductivity().unit == u.W / u.m / u.K)
        
    def test_ion_viscosity_units(self):
        """output should be Quantity with units of Pa s """
        assert(self.ct.ion_viscosity().unit == u.Pa * u.s)
    
    def test_electron_viscosity_units(self):
        """output should be Quantity with units of Pa s"""
        assert(self.ct.electron_viscosity().unit == u.Pa * u.s)

    
#for Braginskii model, from Ref 1:
#nondim_coefficients
#for hall >> 1, (100 or 1000 should do it):
#
#from Table 1:
#eq (2.8), Z = 1, alpha_e_par_hat = 0.51
#eq (2.8), Z = 2, alpha_e_par_hat = 0.44
#eq (2.8), Z = 3, alpha_e_par_hat = 0.40
#eq (2.8), Z = 4, alpha_e_par_hat = 0.38
#eq (2.8), Z = inf, alpha_e_par_hat = 0.29
#eq (2.9), Z = 1, beta_e_par_hat = 0.71
#eq (2.9), Z = 2, beta_e_par_hat = 0.9
#eq (2.9), Z = 3, beta_e_par_hat = 1.0
#eq (2.9), Z = 4, beta_e_par_hat = 1.1
#eq (2.9), Z = inf, beta_e_par_hat = 1.5
#eq (2.12), Z = 1, kappa_e_par_hat = 3.16
#eq (2.12), Z = 2, kappa_e_par_hat = 4.9
#eq (2.12), Z = 3, kappa_e_par_hat = 6.1
#eq (2.12), Z = 4, kappa_e_par_hat = 6.9
#eq (2.12), Z = inf, kappa_e_par_hat = 12.5
#eq (2.13), Z = 1, kappa_e_perp_hat * hall^2 = 4.66
#eq (2.13), Z = 2, kappa_e_perp_hat * hall^2 = 4.0
#eq (2.13), Z = 3, kappa_e_perp_hat * hall^2 = 3.7
#eq (2.13), Z = 4, kappa_e_perp_hat * hall^2 = 3.6
#eq (2.13), Z = inf, kappa_e_perp_hat * hall^2 = 3.2
#
#eq (2.15), for all Z, kappa_i_par_hat = 3.9
#eq (2.16), for all Z, kappa_i_perp_hat * hall^2 = 2.0
#eq (2.22), eta_i_hat_0 = 0.96
#eq (2.23), eta_i_hat_1 * hall^2 = 0.3
#eq (2.23), eta_i_hat_2 * hall^2 = 1.2
#eq (2.24), eta_i_hat_3 * hall = 0.5
#eq (2.24), eta_i_hat_4 * hall = 1.0
#eq (2.25), eta_e_hat_0 = 0.73
#eq (2.26), eta_e_hat_1 * hall^2 = 0.51
#eq (2.26), eta_e_hat_2 * hall^2 = 2.04
#eq (2.27), eta_e_hat_3 * hall = 0.5
#eq (2.27), eta_e_hat_4 * hall = 1.0


# test class for _nondim_tc_e_braginskii function:
class Test__nondim_tc_e_braginskii(object):
    
    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
            
    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'par', 3.16), # eq (2.12), table 1
        (2, 'par', 4.9),  # eq (2.12), table 1
        (3, 'par', 6.1),  # eq (2.12), table 1
        (4, 'par', 6.9),  # eq (2.12), table 1
        (np.inf, 'par', 12.5), # eq (2.12), table 1
    ])
    def test_known_values_par(self, Z, field_orientation, expected):
        kappa_e_hat = _nondim_tc_e_braginskii(self.big_hall, 
                                              Z, 
                                              field_orientation) 
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_e_hat, decimal_places) == expected)
        
    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'perp', 4.66), # eq (2.13), table 1
        (2, 'perp', 4.0),  # eq (2.13), table 1
        (3, 'perp', 3.7),  # eq (2.13), table 1
        (4, 'perp', 3.6),  # eq (2.13), table 1
        (np.inf, 'perp', 3.2),  # eq (2.13),table 1
    ])
    def test_known_values_perp(self, Z, field_orientation, expected):
        kappa_e_hat = _nondim_tc_e_braginskii(self.big_hall, 
                                              Z, 
                                              field_orientation)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_e_hat * self.big_hall ** 2, 
                        decimal_places) == expected)


# test class for _nondim_tc_i_braginskii function:
class Test__nondim_tc_i_braginskii(object):
    
    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
            
    def test_known_values_par(self):
        kappa_i_hat = _nondim_tc_i_braginskii(self.big_hall, 
                                              field_orientation='par') 
        expected = 3.9 # Braginskii '65 eq (2.15)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_i_hat, decimal_places) == expected)
        
    def test_known_values_perp(self):
        kappa_i_hat = _nondim_tc_i_braginskii(self.big_hall, 
                                              field_orientation='perp')
        expected = 2.0 # Braginskii '65 eq (2.16)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_i_hat * self.big_hall ** 2, 
                        decimal_places) == expected)


# test class for _nondim_tec_braginskii function:
class Test__nondim_tec_braginskii(object):
    
    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        
    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'par', 0.71), # eq (2.9), table 1
        (2, 'par', 0.9),  # eq (2.9), table 1
        (3, 'par', 1.0),  # eq (2.9), table 1
        (4, 'par', 1.1),  # eq (2.9), table 1
        (np.inf, 'par', 1.5),  # eq (2.9),table 1
    ])
    def test_known_values_par(self, Z, field_orientation, expected):
        beta_hat = _nondim_tec_braginskii(self.big_hall, 
                                              Z, 
                                              field_orientation) 
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(beta_hat, decimal_places) == expected)


# test class for _nondim_resist_braginskii function:
class Test__nondim_resist_braginskii(object):
    
    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        
    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'par', 0.51), # eq (2.8), table 1
        pytest.param(2, 'par', 0.44, marks=pytest.mark.xfail),  # table 1
        (3, 'par', 0.40),  # eq (2.8), table 1
        (4, 'par', 0.38),  # eq (2.8), table 1
        (np.inf, 'par', 0.29),  # eq (2.8),table 1
    ])
    def test_known_values_par(self, Z, field_orientation, expected):
        beta_hat = _nondim_resist_braginskii(self.big_hall, 
                                              Z, 
                                              field_orientation) 
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(beta_hat, decimal_places) == expected)
        
        
# test class for _nondim_visc_i_braginskii function:
class Test__nondim_visc_i_braginskii(object):
    
    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        
    # values from Braginskii '65
    @pytest.mark.parametrize("expected, idx", [
        (0.96, 0), # eq (2.22)
        (0.3, 1),  # eq (2.23)
        (1.2, 2),  # eq (2.23)
        (0.5, 3),  # eq (2.24)
        (1.0, 4),  # eq (2.24)
    ])
    def test_known_values(self, expected, idx):
        beta_hat = _nondim_visc_i_braginskii(self.big_hall) 
        decimal_places = count_decimal_places(str(expected))
        if idx == 0:
            assert(np.round(beta_hat[idx], decimal_places) == expected)
        elif idx == 1 or idx == 2:
            assert(np.round(beta_hat[idx] * self.big_hall ** 2,
                            decimal_places) == expected)
        elif idx == 3 or idx == 4:
            assert(np.round(beta_hat[idx] * self.big_hall,
                            decimal_places) == expected)


# test class for _nondim_visc_e_braginskii function:
class Test__nondim_visc_e_braginskii(object):
    
    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        
    # values from Braginskii '65
    @pytest.mark.parametrize("Z, expected, idx", [
        (1, 0.73, 0), # eq (2.25)
        (1, 0.51, 1), # eq (2.26)
        pytest.param(1, 2.04, 2, marks=pytest.mark.xfail), # eq (2.26)
        (1, 0.5, 3), # eq (2.27)
        (1, 1.0, 4), # eq (2.27)
    ])
    def test_known_values(self, Z, expected, idx):
        beta_hat = _nondim_visc_e_braginskii(self.big_hall, Z) 
        decimal_places = count_decimal_places(str(expected))
        if idx == 0:
            assert(np.round(beta_hat[idx], decimal_places) == expected)
        elif idx == 1 or idx == 2:
            assert(np.round(beta_hat[idx] * self.big_hall ** 2,
                            decimal_places) == expected)
        elif idx == 3 or idx == 4:
            assert(np.round(beta_hat[idx] * self.big_hall,
                            decimal_places) == expected)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
