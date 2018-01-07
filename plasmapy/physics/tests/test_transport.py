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
                         _nondim_tc_e_spitzer,
                         _nondim_tec_spitzer,
                         _nondim_resist_spitzer,
                         _nondim_tc_e_ji_held,
                         _nondim_tc_i_ji_held,
                         _nondim_tec_ji_held,
                         _nondim_resist_ji_held,
                         _nondim_visc_e_ji_held,
                         _nondim_visc_e_ji_held,
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
            T_e=self.T_e,
            n_e=self.n_e,
            T_i=self.T_i,
            n_i=self.n_i,
            ion_particle=self.ion_particle,
            Z=self.Z,
            B=self.B,
            model=self.model,
            field_orientation=self.field_orientation,
            coulomb_log_ei=self.coulomb_log_val_ei,
            coulomb_log_ii=self.coulomb_log_val_ii,
            V_ei=self.V_ei,
            V_ii=self.V_ii,
            hall_e=self.hall_e,
            hall_i=self.hall_i,
            mu=self.mu,
            theta=self.theta,
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


# test class for _nondim_tc_e_braginskii function:
class Test__nondim_tc_e_braginskii(object):

    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        self.small_hall = 0

    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'par', 3.16),  # eq (2.12), table 1
        (2, 'par', 4.9),  # eq (2.12), table 1
        (3, 'par', 6.1),  # eq (2.12), table 1
        (4, 'par', 6.9),  # eq (2.12), table 1
        (np.inf, 'par', 12.5),  # eq (2.12), table 1
    ])
    def test_known_values_par(self, Z, field_orientation, expected):
        """check some known values"""
        kappa_e_hat = _nondim_tc_e_braginskii(self.big_hall,
                                              Z,
                                              field_orientation)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_e_hat, decimal_places) == expected)

    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'perp', 4.66),  # eq (2.13), table 1
        (2, 'perp', 4.0),  # eq (2.13), table 1
        (3, 'perp', 3.7),  # eq (2.13), table 1
        (4, 'perp', 3.6),  # eq (2.13), table 1
        (np.inf, 'perp', 3.2),  # eq (2.13),table 1
    ])
    def test_known_values_perp(self, Z, field_orientation, expected):
        """check some known values"""
        kappa_e_hat = _nondim_tc_e_braginskii(self.big_hall,
                                              Z,
                                              field_orientation)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_e_hat * self.big_hall ** 2,
                        decimal_places) == expected)

    @pytest.mark.parametrize("Z", [1, 2, 3, 4, np.inf])    
    def test_unmagnetized(self, Z):
        """confirm perp -> par as B -> 0"""
        kappa_e_hat_par = _nondim_tc_e_braginskii(self.small_hall, Z, 'par')
        kappa_e_hat_perp = _nondim_tc_e_braginskii(self.small_hall, Z, 'perp')
        assert np.isclose(kappa_e_hat_par, kappa_e_hat_perp, rtol=1e-3)


# test class for _nondim_tc_i_braginskii function:
class Test__nondim_tc_i_braginskii(object):

    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        self.small_hall = 0

    def test_known_values_par(self):
        """check some known values"""
        kappa_i_hat = _nondim_tc_i_braginskii(self.big_hall,
                                              field_orientation='par')
        expected = 3.9  # Braginskii '65 eq (2.15)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_i_hat, decimal_places) == expected)

    def test_known_values_perp(self):
        """check some known values"""
        kappa_i_hat = _nondim_tc_i_braginskii(self.big_hall,
                                              field_orientation='perp')
        expected = 2.0  # Braginskii '65 eq (2.16)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(kappa_i_hat * self.big_hall ** 2,
                        decimal_places) == expected)
        
    def test_unmagnetized(self):
        """confirm perp -> par as B -> 0"""
        kappa_i_hat_par = _nondim_tc_i_braginskii(self.small_hall, 'par')
        kappa_i_hat_perp = _nondim_tc_i_braginskii(self.small_hall, 'perp')
        assert np.isclose(kappa_i_hat_par, kappa_i_hat_perp, rtol=1e-3)


# test class for _nondim_tec_braginskii function:
class Test__nondim_tec_braginskii(object):

    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        self.small_hall = 0

    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'par', 0.71),  # eq (2.9), table 1
        (2, 'par', 0.9),  # eq (2.9), table 1
        (3, 'par', 1.0),  # eq (2.9), table 1
        (4, 'par', 1.1),  # eq (2.9), table 1
        (np.inf, 'par', 1.5),  # eq (2.9),table 1
    ])
    def test_known_values_par(self, Z, field_orientation, expected):
        """check some known values"""
        beta_hat = _nondim_tec_braginskii(self.big_hall,
                                          Z,
                                          field_orientation)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(beta_hat, decimal_places) == expected)
    
    @pytest.mark.parametrize("Z", [1, 2, 3, 4, np.inf])  
    def test_unmagnetized(self, Z):
        """confirm perp -> par as B -> 0"""
        beta_hat_par = _nondim_tec_braginskii(self.small_hall, Z, 'par')
        beta_hat_perp = _nondim_tec_braginskii(self.small_hall, Z, 'perp')
        assert np.isclose(beta_hat_par, beta_hat_perp, rtol=1e-3)


# test class for _nondim_resist_braginskii function:
class Test__nondim_resist_braginskii(object):

    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        self.small_hall = 0

    # values from Braginskii '65
    @pytest.mark.parametrize("Z, field_orientation, expected", [
        (1, 'par', 0.51),  # eq (2.8), table 1
        pytest.param(2, 'par', 0.44, marks=pytest.mark.xfail),  # table 1
        (3, 'par', 0.40),  # eq (2.8), table 1
        (4, 'par', 0.38),  # eq (2.8), table 1
        (np.inf, 'par', 0.29),  # eq (2.8),table 1
    ])
    def test_known_values_par(self, Z, field_orientation, expected):
        """check some known values"""
        beta_hat = _nondim_resist_braginskii(self.big_hall,
                                             Z,
                                             field_orientation)
        decimal_places = count_decimal_places(str(expected))
        assert(np.round(beta_hat, decimal_places) == expected)
    
    @pytest.mark.parametrize("Z", [1, 2, 3, 4, np.inf])  
    def test_unmagnetized(self, Z):
        """confirm perp -> par as B -> 0"""
        alpha_hat_par = _nondim_resist_braginskii(self.small_hall, Z, 'par')
        alpha_hat_perp = _nondim_resist_braginskii(self.small_hall, Z, 'perp')
        assert np.isclose(alpha_hat_par, alpha_hat_perp, rtol=1e-3)


# test class for _nondim_visc_i_braginskii function:
class Test__nondim_visc_i_braginskii(object):

    def setup_method(self):
        """set up some initial values for tests"""
        self.big_hall = 1000
        self.small_hall = 0

    # values from Braginskii '65
    @pytest.mark.parametrize("expected, idx", [
        (0.96, 0),  # eq (2.22)
        (0.3, 1),  # eq (2.23)
        (1.2, 2),  # eq (2.23)
        (0.5, 3),  # eq (2.24)
        (1.0, 4),  # eq (2.24)
    ])
    def test_known_values(self, expected, idx):
        """check some known values"""
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
        self.small_hall = 0

    # values from Braginskii '65
    @pytest.mark.parametrize("Z, expected, idx", [
        (1, 0.73, 0),  # eq (2.25)
        (1, 0.51, 1),  # eq (2.26)
        pytest.param(1, 2.04, 2, marks=pytest.mark.xfail),  # eq (2.26)
        (1, 0.5, 3),  # eq (2.27)
        (1, 1.0, 4),  # eq (2.27)
    ])
    def test_known_values(self, Z, expected, idx):
        """check some known values"""
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


@pytest.mark.parametrize("Z", [1, 2, 4, 16, np.inf])
def test__nondim_tc_e_spitzer(Z):
    """test _nondim_tc_e_spitzer function"""
    kappa = _nondim_tc_e_spitzer(Z)
    if Z == 1:
        kappa_check = 3.203
        rtol = 1e-3
    elif Z == 2 or Z == 4:
        kappa_check = _nondim_tc_e_braginskii(0, Z, 'par')
        rtol = 2e-2
    elif Z == 16:
        kappa_check = _nondim_tc_e_ji_held(0, Z, 'par')
        rtol = 2e-2
    elif Z == np.inf:
        kappa_check = _nondim_tc_e_ji_held(0, 1e6, 'par')
        rtol = 2e-2
    assert np.isclose(kappa, kappa_check, rtol=rtol)


@pytest.mark.parametrize("Z", [1, 2, 4, 16, np.inf])
def test__nondim_resist_spitzer(Z):
    """test _nondim_resist_spitzer function"""
    alpha = _nondim_resist_spitzer(Z)
    if Z == 1:
        alpha_check = 0.5064
        rtol = 1e-3
    elif Z == 2 or Z == 4 or Z == np.inf:
        alpha_check = _nondim_resist_braginskii(0, Z, 'par')
        rtol = 2e-2
    elif Z == 16:
        alpha_check = _nondim_resist_ji_held(0, Z, 'par')
        rtol = 2e-2
    assert np.isclose(alpha, alpha_check, rtol=rtol)


@pytest.mark.parametrize("Z", [1, 2, 4, 16, np.inf]) 
def test__nondim_tec_spitzer(Z):
    """test _nondim_tec_spitzer function"""
    beta = _nondim_tec_spitzer(Z)
    if Z == 1:
        beta_check = 0.699
        rtol = 1e-3
    elif Z == 2 or Z == 4 or Z == np.inf:
        beta_check = _nondim_tec_braginskii(0, Z, 'par')
        rtol = 2e-2
    elif Z == 16:
        beta_check = _nondim_tec_ji_held(0, Z, 'par')
        rtol = 2e-2
    assert np.isclose(beta, beta_check, rtol=rtol)