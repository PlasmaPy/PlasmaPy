import astropy.units as u
import numpy as np
import pytest

from astropy.constants import c

from plasmapy.formulary.quantum import (
    _chemical_potential_interp,
    chemical_potential,
    deBroglie_wavelength,
    Ef_,
    Fermi_energy,
    lambdaDB_,
    lambdaDB_th_,
    quantum_theta,
    thermal_deBroglie_wavelength,
    Thomas_Fermi_length,
    Wigner_Seitz_radius,
)
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import RelativityError


def test_deBroglie_wavelength():

    dbwavelength1 = deBroglie_wavelength(2e7 * u.cm / u.s, "e")
    assert np.isclose(dbwavelength1.value, 3.628845222852886e-11)
    assert dbwavelength1.unit == u.m

    dbwavelength2 = deBroglie_wavelength(0 * u.m / u.s, "e")
    assert dbwavelength2 == np.inf * u.m

    V_array = np.array([2e5, 0]) * u.m / u.s
    dbwavelength_arr = deBroglie_wavelength(V_array, "e")

    assert np.isclose(dbwavelength_arr.value[0], 3.628845222852886e-11)
    assert dbwavelength_arr.value[1] == np.inf
    assert dbwavelength_arr.unit == u.m

    V_array = np.array([2e5, 2e5]) * u.m / u.s
    dbwavelength_arr = deBroglie_wavelength(V_array, "e")

    assert np.isclose(dbwavelength_arr.value[0], 3.628845222852886e-11)
    assert np.isclose(dbwavelength_arr.value[1], 3.628845222852886e-11)
    assert dbwavelength_arr.unit == u.m

    assert deBroglie_wavelength(-5e5 * u.m / u.s, "p") == deBroglie_wavelength(
        5e5 * u.m / u.s, "p"
    )

    assert deBroglie_wavelength(-5e5 * u.m / u.s, "e+") == deBroglie_wavelength(
        5e5 * u.m / u.s, "e"
    )

    assert deBroglie_wavelength(1 * u.m / u.s, 5 * u.kg) == deBroglie_wavelength(
        100 * u.cm / u.s, 5000 * u.g
    )


@pytest.mark.parametrize(
    "kwargs, exception",
    [
        ({"V": c * 1.00000001, "particle": "e-"}, RelativityError),
        ({"V": 8 * u.m / u.s, "particle": 5 * u.m}, u.UnitConversionError),
        ({"V": 8 * u.m / u.s, "particle": "invalid particle"}, InvalidParticleError),
    ],
)
def test_deBroglie_exceptions(kwargs, exception):
    with pytest.raises(exception):
        deBroglie_wavelength(**kwargs)


def test_deBroglie_warning():
    with pytest.warns(u.UnitsWarning):
        deBroglie_wavelength(0.79450719277, "Be-7 1+")


# defining some plasma parameters for tests
T_e = 1 * u.eV
n_e = 1e23 * u.cm**-3
# should probably change this to use unittest module
# add tests for numpy arrays as inputs
# add tests for different astropy units (random fuzzing method?)


def test_thermal_deBroglie_wavelength():
    r"""Test the thermal_deBroglie_wavelength function in quantum.py."""
    lambda_dbTh = thermal_deBroglie_wavelength(T_e)
    # true value at 1 eV
    lambda_dbTh_true = 6.919367518364532e-10
    # test a simple case for expected value
    expectStr = (
        "Thermal deBroglie wavelength at 1 eV should be "
        f"{lambda_dbTh_true} and not {lambda_dbTh}"
    )
    assert np.isclose(
        lambda_dbTh.value, lambda_dbTh_true, rtol=1e-5, atol=0.0
    ), expectStr
    # testing returned units
    assert lambda_dbTh.unit == u.m
    # testing exceptions
    with pytest.raises(TypeError):
        thermal_deBroglie_wavelength("Bad Input")
    with pytest.raises(ValueError):
        thermal_deBroglie_wavelength(T_e=-1 * u.eV)


def test_Fermi_energy():
    r"""Test the Fermi_energy function in quantum.py."""
    energy_F = Fermi_energy(n_e)
    # true value at 1e23 cm-3
    energy_F_true = 1.2586761116196002e-18
    # test a simple case for expected value
    expectStr = (
        "Fermi energy at 1e23 cm^-3 should be {energy_F_true} and not {energy_F}."
    )
    assert np.isclose(energy_F.value, energy_F_true, rtol=1e-5, atol=0.0), expectStr
    # testing returned units
    assert energy_F.unit == u.J
    # testing exceptions
    with pytest.raises(TypeError):
        Fermi_energy("Bad Input")
    with pytest.raises(ValueError):
        Fermi_energy(n_e=-1 * u.m**-3)


def test_Thomas_Fermi_length():
    r"""Test the Thomas_Fermi_length function in quantum.py."""
    lambda_TF = Thomas_Fermi_length(n_e)
    # true value at 1e23 cm-3
    lambda_TF_true = 5.379914085596706e-11
    # test a simple case for expected value
    expectStr = (
        "Thomas-Fermi length at 1e23 cm^-3 should be "
        f"{lambda_TF_true} and not {lambda_TF}."
    )
    assert np.isclose(lambda_TF.value, lambda_TF_true, rtol=1e-5, atol=0.0), expectStr
    # testing returned units
    assert lambda_TF.unit == u.m
    # testing exceptions
    with pytest.raises(TypeError):
        Thomas_Fermi_length("Bad Input")
    with pytest.raises(ValueError):
        Thomas_Fermi_length(n_e=-1 * u.m**-3)


def test_Wigner_Seitz_radius():
    """
    Checks Wigner-Seitz radius for a known value.
    """
    n_e = 1e23 * u.cm**-3
    radiusTrue = 1.3365046175719772e-10 * u.m
    radiusMeth = Wigner_Seitz_radius(n_e)
    testTrue = u.isclose(radiusMeth, radiusTrue, rtol=1e-5)
    errStr = f"Error in Wigner_Seitz_radius(), got {radiusMeth}, should be {radiusTrue}"
    assert testTrue, errStr


@pytest.mark.slow
class TestChemicalPotential:

    value_test_parameters = (
        "n_e, T, expected_value",
        [
            (1e20 * u.cm**-3, 11604 * u.K, 0),
            (1e25 * u.cm**-3, 11604 * u.K, 268.68166791746324),
        ],
    )

    @staticmethod
    @pytest.mark.parametrize(*value_test_parameters)
    def test_return_value(n_e, T, expected_value):
        """
        Tests chemical_potential for expected value.
        """
        calculated_value = chemical_potential(n_e, T)

        error_message = f"Chemical potential value should be {expected_value} and not {calculated_value.value}."
        assert np.isclose(
            calculated_value.value, expected_value, rtol=1e-16, atol=0.0
        ), error_message
        assert calculated_value.unit == u.dimensionless_unscaled

    @staticmethod
    @pytest.mark.parametrize(*value_test_parameters)
    def test_fail1(n_e, T, expected_value):
        """
        Tests if test_return_value would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        expected_failure_value = expected_value + 1e-10
        calculated_value = chemical_potential(n_e, T)

        error_message = (
            f"Chemical potential value test gives {calculated_value.value} and "
            f"should not be equal to {expected_failure_value}."
        )
        assert not u.isclose(
            calculated_value.value, expected_failure_value, rtol=1e-16, atol=0.0
        ), error_message


class Test__chemical_potential_interp:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.n_e = 1e23 * u.cm**-3
        cls.T = 11604 * u.K
        cls.True1 = 7.741254037813922

    def test_known1(self):
        """
        Tests Fermi_integral for expected value.
        """
        methodVal = _chemical_potential_interp(self.n_e, self.T)
        testTrue = u.isclose(methodVal, self.True1, rtol=1e-16, atol=0.0)
        errStr = f"Chemical potential value should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 + 1e-15
        methodVal = _chemical_potential_interp(self.n_e, self.T)
        testTrue = not u.isclose(methodVal, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Chemical potential value test gives {methodVal} "
            f"and should not be equal to {fail1}."
        )
        assert testTrue, errStr


def test_quantum_aliases():
    r"""Test all aliases defined in quantum.py"""

    assert Ef_ is Fermi_energy
    assert lambdaDB_ is deBroglie_wavelength
    assert lambdaDB_th_ is thermal_deBroglie_wavelength


class TestQuantumTheta:
    """Test the quantum_theta function in quantum.py."""

    def test_units(self):
        """Test the return units"""

        theta = quantum_theta(1 * u.eV, 1e26 * u.m**-3)

        assert theta.unit.is_equivalent(u.dimensionless_unscaled)

    @pytest.mark.parametrize(
        "T, n_e, expected_theta",
        [
            (1 * u.eV, 1e26 * u.m**-3, 12.72906),  # Both regimes are present
            (2 / 3 * u.eV, 1e40 * u.m**-3, 3.93887e-9),  # Fermi regime
            (8.61e2 * u.eV, 1e12 * u.m**-3, 2.36120e13),  # Thermal regime
            (1 * u.K, 1e26 * u.m**-3, 1.09690e-3),  # Specify temperature in Kelvin
        ],
    )
    def test_value(self, T, n_e, expected_theta):
        """Compare the calculated theta with the expected value."""

        theta = quantum_theta(T, n_e)

        assert np.isclose(theta.value, expected_theta)
