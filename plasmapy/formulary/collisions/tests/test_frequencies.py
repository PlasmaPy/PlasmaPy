import astropy.units as u
import numpy as np
import pytest

from astropy.constants import k_B, m_p

from plasmapy.formulary.collisions.frequencies import (
    collision_frequency,
    fundamental_electron_collision_freq,
    fundamental_ion_collision_freq,
    MaxwellianCollisionFrequencies,
    SingleParticleCollisionFrequencies,
)
from plasmapy.particles import Particle
from plasmapy.utils import exceptions
from plasmapy.utils.exceptions import CouplingWarning, PhysicsError
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray


class TestSingleParticleCollisionFrequencies:
    """Test the SingleParticleCollisionFrequencies class in collisions.py."""

    attribute_units_test_case = SingleParticleCollisionFrequencies(
        Particle("e-"),
        Particle("e-"),
        v_drift=1 * u.m / u.s,
        T_b=1 * u.K,
        n_b=1 * u.m**-3,
        Coulomb_log=1,
    )

    MKS_unit_conversion_test_constructor_arguments = {
        "test_particle": Particle("e-"),
        "field_particle": Particle("e-"),
        "v_drift": 1e5 * u.m / u.s,
        "T_b": 1e3 * u.eV,
        "n_b": 1e26 * u.m**-3,
        "Coulomb_log": 10 * u.dimensionless_unscaled,
    }

    arguments_to_convert = ["v_drift", "n_b"]

    CGS_unit_conversion_test_constructor_arguments = (
        MKS_unit_conversion_test_constructor_arguments
    )

    for argument_to_convert in arguments_to_convert:
        CGS_unit_conversion_test_constructor_arguments[
            argument_to_convert
        ] = CGS_unit_conversion_test_constructor_arguments[argument_to_convert].cgs

    MKS_test_case = SingleParticleCollisionFrequencies(
        **MKS_unit_conversion_test_constructor_arguments
    )
    CGS_test_case = SingleParticleCollisionFrequencies(
        **CGS_unit_conversion_test_constructor_arguments
    )

    return_values_to_test = [
        "momentum_loss",
        "transverse_diffusion",
        "parallel_diffusion",
        "energy_loss",
    ]

    ones_array = np.ones(5)
    ones_array2d = np.ones([5, 5])

    @pytest.mark.parametrize(
        "attribute_to_test, expected_attribute_units",
        [
            ("momentum_loss", u.Hz),
            ("transverse_diffusion", u.Hz),
            ("parallel_diffusion", u.Hz),
            ("energy_loss", u.Hz),
            ("x", u.dimensionless_unscaled),
            ("Lorentz_collision_frequency", u.Hz),
            ("Coulomb_log", u.dimensionless_unscaled),
        ],
    )
    def test_units(self, attribute_to_test, expected_attribute_units):
        """Test the return units"""

        assert getattr(
            self.attribute_units_test_case, attribute_to_test
        ).unit.is_equivalent(expected_attribute_units)

    @pytest.mark.parametrize(
        "attribute_to_test",
        [
            "momentum_loss",
            "transverse_diffusion",
            "parallel_diffusion",
            "energy_loss",
            "x",
            "Lorentz_collision_frequency",
        ],
    )
    def test_conversion_consistency(self, attribute_to_test):
        """
        Test that a consistent value is computed for attributes
        regardless of argument units.
        """

        MKS_result = getattr(self.MKS_test_case, attribute_to_test)
        CGS_result = getattr(self.CGS_test_case, attribute_to_test)

        assert MKS_result == CGS_result

    @staticmethod
    def get_limit_value(interaction_type, limit_type, cases):
        """
        Get the limiting values for frequencies given the two particles
        interacting, and their frequencies class.

        These formulae are taken from page 31 of the NRL Formulary.
        """

        v_a = (0.5 * cases.test_particle.mass * cases.v_drift**2).to(u.eV).value
        T_b = (cases.T_b * k_B).to(u.eV).value

        limit_values = []

        if interaction_type == "e|e":
            if limit_type == "slow":
                limit_values.extend(
                    [
                        5.8e-6 * T_b ** (-1.5),
                        5.8e-6 * T_b ** (-0.5) * v_a ** (-1),
                        2.9e-6 * T_b ** (-0.5) * v_a ** (-1),
                    ]
                )
            elif limit_type == "fast":
                limit_values.extend(
                    [
                        7.7e-6 * v_a ** (-1.5),
                        7.7e-6 * v_a ** (-1.5),
                        3.9e-6 * T_b * v_a ** (-2.5),
                    ]
                )
        elif interaction_type == "e|i":
            mu = (cases.field_particle.mass / m_p).value

            if limit_type == "slow":
                limit_values.extend(
                    [
                        0.23 * mu**1.5 * T_b**-1.5,
                        2.5e-4 * mu**0.5 * T_b**-0.5 * v_a**-1,
                        1.2e-4 * mu**0.5 * T_b**-0.5 * v_a**-1,
                    ]
                )
            elif limit_type == "fast":
                limit_values.extend(
                    [
                        3.9e-6 * v_a**-1.5,
                        7.7e-6 * v_a**-1.5,
                        2.1e-9 * mu**-1 * T_b * v_a**-2.5,
                    ]
                )
        elif interaction_type == "i|e":
            mu = (cases.test_particle.mass / m_p).value

            if limit_type == "slow":
                limit_values.extend(
                    [
                        1.6e-9 * mu**-1 * T_b ** (-1.5),
                        3.2e-9 * mu**-1 * T_b ** (-0.5) * v_a**-1,
                        1.6e-9 * mu**-1 * T_b ** (-0.5) * v_a**-1,
                    ]
                )
            elif limit_type == "fast":
                limit_values.extend(
                    [
                        1.7e-4 * mu**0.5 * v_a**-1.5,
                        1.8e-7 * mu**-0.5 * v_a**-1.5,
                        1.7e-4 * mu**0.5 * T_b * v_a**-2.5,
                    ]
                )

        elif interaction_type == "i|i":
            mu = (cases.test_particle.mass / m_p).value
            mu_prime = (cases.field_particle.mass / m_p).value

            if limit_type == "slow":
                limit_values.extend(
                    [
                        6.8e-8
                        * mu_prime**0.5
                        * mu**-1
                        * (1 + mu_prime / mu)
                        * T_b**-1.5,
                        1.4e-7 * mu_prime**0.5 * mu**-1 * T_b**-0.5 * v_a**-1,
                        6.8e-8 * mu_prime**0.5 * mu**-1 * T_b**-0.5 * v_a**-1,
                    ]
                )
            elif limit_type == "fast":
                limit_values.extend(
                    [
                        9e-8 * (1 / mu + 1 / mu_prime) * mu**0.5 * v_a**-1.5,
                        1.8 * 10**-7 * mu**-0.5 * v_a**-1.5,
                        9e-8 * mu**0.5 * mu_prime**-1 * T_b * v_a**-2.5,
                    ]
                )
        # The expected energy loss collision frequency should always equal this
        limit_values.append(
            2 * cases.momentum_loss.value
            - cases.transverse_diffusion.value
            - cases.parallel_diffusion.value
        )

        return limit_values

    @pytest.mark.parametrize(
        "interaction_type, limit_type, constructor_arguments, "
        "constructor_keyword_arguments",
        [
            # Slow limit (x << 1)
            (
                "e|e",
                "slow",
                (Particle("e-"), Particle("e-")),
                {
                    "v_drift": 1 * u.cm / u.s,
                    "T_b": 1e4 * u.eV,
                    "n_b": 1e15 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "e|i",
                "slow",
                (Particle("e-"), Particle("Na+")),
                {
                    "v_drift": 1 * u.cm / u.s,
                    "T_b": 1e4 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "e|i",
                "slow",
                (Particle("e-"), Particle("Ba 2+")),
                {
                    "v_drift": 1 * u.cm / u.s,
                    "T_b": 1e4 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|e",
                "slow",
                (Particle("Na+"), Particle("e-")),
                {
                    "v_drift": 1 * u.cm / u.s,
                    "T_b": 1e2 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|e",
                "slow",
                (Particle("Be 2+"), Particle("e-")),
                {
                    "v_drift": 1 * u.cm / u.s,
                    "T_b": 1e2 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|i",
                "slow",
                (Particle("Na+"), Particle("Cl-")),
                {
                    "v_drift": 1 * u.cm / u.s,
                    "T_b": 1e4 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|i",
                "slow",
                (Particle("Na+"), Particle("S 2-")),
                {
                    "v_drift": 1 * u.cm / u.s,
                    "T_b": 1e4 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            # Fast limit (x >> 1)
            (
                "e|e",
                "fast",
                (Particle("e-"), Particle("e-")),
                {
                    "v_drift": 6e8 * u.cm / u.s,
                    "T_b": 1e-1 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "e|i",
                "fast",
                (Particle("e-"), Particle("Na+")),
                {
                    "v_drift": 6e5 * u.cm / u.s,
                    "T_b": 1e-3 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "e|i",
                "fast",
                (Particle("e-"), Particle("Zn 2+")),
                {
                    "v_drift": 6e5 * u.cm / u.s,
                    "T_b": 1e-3 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|e",
                "fast",
                (Particle("Na+"), Particle("e-")),
                {
                    "v_drift": 3e7 * u.cm / u.s,
                    "T_b": 1e-3 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|e",
                "fast",
                (Particle("Ca 2+"), Particle("e-")),
                {
                    "v_drift": 3e7 * u.cm / u.s,
                    "T_b": 1e-3 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|i",
                "fast",
                (Particle("Na+"), Particle("Cl-")),
                {
                    "v_drift": 3e7 * u.cm / u.s,
                    "T_b": 1e2 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "i|i",
                "fast",
                (Particle("Be 2+"), Particle("Cl-")),
                {
                    "v_drift": 3e7 * u.cm / u.s,
                    "T_b": 1e2 * u.eV,
                    "n_b": 1e20 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
        ],
    )
    def test_limit_values(
        self,
        interaction_type,
        limit_type,
        constructor_arguments,
        constructor_keyword_arguments,
    ):
        """Test the return values"""

        value_test_case = SingleParticleCollisionFrequencies(
            *constructor_arguments, **constructor_keyword_arguments
        )

        coulomb_density_constant = (
            constructor_keyword_arguments["Coulomb_log"].value
            * constructor_keyword_arguments["n_b"].to(u.cm**-3).value
        )

        expected_limit_values = self.get_limit_value(
            interaction_type, limit_type, value_test_case
        )

        if interaction_type == "e|e":
            charge_constant = 1
        elif interaction_type == "e|i":
            charge_constant = value_test_case.field_particle.charge_number**2
        elif interaction_type == "i|e":
            charge_constant = value_test_case.test_particle.charge_number**2
        elif interaction_type == "i|i":
            charge_constant = (
                value_test_case.test_particle.charge_number
                * value_test_case.field_particle.charge_number
            ) ** 2

        for attribute_name, expected_limit_value in zip(
            self.return_values_to_test, expected_limit_values
        ):
            calculated_limit_value = getattr(value_test_case, attribute_name).value
            # Energy loss limit value is already in units of
            # frequencies because of the way it is calculated
            if attribute_name != "energy_loss":
                calculated_limit_value = calculated_limit_value / (
                    coulomb_density_constant * charge_constant
                )

            assert np.allclose(
                calculated_limit_value, expected_limit_value, rtol=0.05, atol=0
            )

    @pytest.mark.parametrize(
        "expected_error, constructor_arguments, constructor_keyword_arguments",
        [
            # Arrays of unequal shape error
            (
                ValueError,
                (Particle("e-"), Particle("e-")),
                {
                    "v_drift": np.ndarray([1, 1]) * u.cm / u.s,
                    "T_b": 1 * u.eV,
                    "n_b": ones_array * u.cm**-3,
                    "Coulomb_log": 1 * u.dimensionless_unscaled,
                },
            ),
        ],
    )
    def test_init_errors(
        self, expected_error, constructor_arguments, constructor_keyword_arguments
    ):
        """Test errors raised in the __init__ function body"""

        with pytest.raises(expected_error):
            SingleParticleCollisionFrequencies(
                *constructor_arguments, **constructor_keyword_arguments
            )

    @pytest.mark.parametrize(
        "constructor_keyword_arguments",
        [
            {
                "test_particle": Particle("e-"),
                "field_particle": Particle("e-"),
                "v_drift": ones_array * u.cm / u.s,
                "T_b": ones_array * u.eV,
                "n_b": ones_array * u.cm**-3,
                "Coulomb_log": ones_array * u.dimensionless_unscaled,
            },
            {
                "test_particle": Particle("e-"),
                "field_particle": Particle("e-"),
                "v_drift": ones_array2d * u.m / u.s,
                "T_b": ones_array2d * u.eV,
                "n_b": ones_array2d * u.cm**-3,
                "Coulomb_log": ones_array2d * u.dimensionless_unscaled,
            },
        ],
    )
    def test_handle_ndarrays(self, constructor_keyword_arguments):
        """Test for ability to handle numpy array quantities"""

        SingleParticleCollisionFrequencies(**constructor_keyword_arguments)


class TestMaxwellianCollisionFrequencies:
    ones_array = np.ones(5)

    @staticmethod
    def get_fundamental_frequency(species, n, T_a, Coulomb_log):
        """
        This special case for computing the fundamental frequencies
        comes from page 33 of the NRL Formulary.  The formulary
        provides limiting cases for the
        `Maxwellian_avg_##_collision_freq` family of attributes in the
        case that T_a ~ T_b.
        """

        # Strip the units from these quantities and ensure they are in CGS units
        n = n.to(u.cm**-3).value
        T_a = T_a.to(u.eV).value

        if species.is_electron:
            return (2.9e-6 * n * Coulomb_log * T_a**-1.5) * u.Hz
        elif species.is_ion:
            mu = (species.mass / m_p).value

            return (4.8e-8 * n * Coulomb_log * T_a**-1.5 * mu**-0.5) * u.Hz

    @pytest.mark.parametrize(
        "expected_error, constructor_arguments, constructor_keyword_arguments",
        [
            # Arrays of unequal shape error
            (
                ValueError,
                (Particle("e-"), Particle("e-")),
                {
                    "v_drift": np.array([1, 1]) * u.m / u.s,
                    "T_a": 1 * u.K,
                    "n_a": 1 * u.m**-3,
                    "T_b": 1 * u.K,
                    "n_b": ones_array * u.m**-3,
                    "Coulomb_log": 1 * u.dimensionless_unscaled,
                },
            ),
        ],
    )
    def test_init_errors(
        self, expected_error, constructor_arguments, constructor_keyword_arguments
    ):
        """Test errors raised in the __init__ function body"""

        with pytest.raises(expected_error):
            MaxwellianCollisionFrequencies(
                *constructor_arguments, **constructor_keyword_arguments
            )

    @pytest.mark.parametrize(
        "expected_error, constructor_arguments, "
        "constructor_keyword_arguments, attribute_name",
        [
            # Specified interaction isn't electron-ion
            (
                ValueError,
                ("e-", "e-"),
                {
                    "v_drift": 1e-5 * u.cm / u.s,
                    "T_a": 1e3 * u.eV,
                    "n_a": 1e10 * u.cm**-3,
                    "T_b": 1e3 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
                "Maxwellian_avg_ei_collision_freq",
            ),
            # Specified interaction isn't ion-ion
            (
                ValueError,
                ("Na+", "e-"),
                {
                    "v_drift": 1e-5 * u.cm / u.s,
                    "T_a": 1e3 * u.eV,
                    "n_a": 1e10 * u.cm**-3,
                    "T_b": 1e3 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
                "Maxwellian_avg_ii_collision_freq",
            ),
            # Populations are not slowly flowing error
            (
                PhysicsError,
                ("e-", "Na+"),
                {
                    "v_drift": 1e10 * u.cm / u.s,
                    "T_a": 1e3 * u.eV,
                    "n_a": 1e10 * u.cm**-3,
                    "T_b": 1e3 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
                "Maxwellian_avg_ei_collision_freq",
            ),
            # Populations are not slowly flowing error
            (
                PhysicsError,
                ("Na+", "Na+"),
                {
                    "v_drift": 1e10 * u.cm / u.s,
                    "T_a": 1e3 * u.eV,
                    "n_a": 1e10 * u.cm**-3,
                    "T_b": 1e3 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
                "Maxwellian_avg_ii_collision_freq",
            ),
        ],
    )
    def test_attribute_errors(
        self,
        expected_error,
        constructor_arguments,
        constructor_keyword_arguments,
        attribute_name,
    ):
        """Test errors raised in attribute bodies"""

        test_case = MaxwellianCollisionFrequencies(
            *constructor_arguments, **constructor_keyword_arguments
        )

        with pytest.raises(expected_error):
            getattr(test_case, attribute_name)

    @pytest.mark.parametrize(
        "frequency_to_test, constructor_keyword_arguments",
        [
            (
                "Maxwellian_avg_ei_collision_freq",
                {
                    "test_particle": Particle("e-"),
                    "field_particle": Particle("Li+"),
                    "v_drift": 1e-5 * u.cm / u.s,
                    "T_a": 1e3 * u.eV,
                    "n_a": 1e10 * u.cm**-3,
                    "T_b": 1e3 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
            (
                "Maxwellian_avg_ii_collision_freq",
                {
                    "test_particle": Particle("Li+"),
                    "field_particle": Particle("Cl-"),
                    "v_drift": 1 * u.cm / u.s,
                    "T_a": 1e3 * u.eV,
                    "n_a": 1e10 * u.cm**-3,
                    "T_b": 1 * u.eV,
                    "n_b": 1e10 * u.cm**-3,
                    "Coulomb_log": 10 * u.dimensionless_unscaled,
                },
            ),
        ],
    )
    def test_fundamental_frequency_values(
        self, frequency_to_test, constructor_keyword_arguments
    ):
        value_test_case = MaxwellianCollisionFrequencies(
            **constructor_keyword_arguments
        )

        calculated_value = getattr(value_test_case, frequency_to_test)
        expected_value = self.get_fundamental_frequency(
            constructor_keyword_arguments["test_particle"],
            constructor_keyword_arguments["n_a"],
            constructor_keyword_arguments["T_a"],
            constructor_keyword_arguments["Coulomb_log"],
        )

        assert np.allclose(calculated_value, expected_value, rtol=5e-3, atol=0)


class Test_collision_frequency:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.n = 1e17 * u.cm**-3
        cls.particles = ("e", "p")
        cls.electrons = ("e", "e")
        cls.protons = ("p", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 1.3468281539854646e12
        cls.True_electrons = 1904702641552.1638
        cls.True_protons = 44450104815.91857
        cls.True_zmean = 1346828153985.4646

    def test_symmetry(self):
        with pytest.warns(CouplingWarning):
            result = collision_frequency(self.T, self.n, self.particles)
            resultRev = collision_frequency(self.T, self.n, self.particles[::-1])
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = collision_frequency(
                self.T,
                self.n,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(self.True1, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = f"Collision frequency should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = collision_frequency(
                self.T,
                self.n,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = not np.isclose(methodVal.si.value, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Collision frequency value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"particles": ("e", "e")},
            {"particles": ("e", "p")},
            {"particles": ("p", "p")},
        ],
    )
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans, kwargs):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            collision_frequency, insert_some_nans, insert_all_nans, kwargs
        )

    def test_electrons(self):
        """
        Testing collision frequency between electrons.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = collision_frequency(
                self.T,
                self.n,
                self.electrons,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(
            self.True_electrons, methodVal.si.value, rtol=1e-1, atol=0.0
        )
        errStr = (
            f"Collision frequency should be {self.True_electrons} and "
            f"not {methodVal}."
        )
        assert testTrue, errStr

    def test_protons(self):
        """
        Testing collision frequency between protons (ions).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = collision_frequency(
                self.T,
                self.n,
                self.protons,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(
            self.True_protons, methodVal.si.value, rtol=1e-1, atol=0.0
        )
        errStr = (
            f"Collision frequency should be {self.True_protons} and "
            f"not {methodVal}."
        )
        assert testTrue, errStr

    def test_zmean(self):
        """
        Test collisional frequency function when given arbitrary z_mean.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = collision_frequency(
                self.T,
                self.n,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(self.True_zmean, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = f"Collision frequency should be {self.True_zmean} and not {methodVal}."
        assert testTrue, errStr


class Test_fundamental_electron_collision_freq:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T_arr = np.array([1, 2]) * u.eV
        cls.n_arr = np.array([1e20, 2e20]) * u.cm**-3
        cls.ion = "p"
        cls.coulomb_log = 10

    # TODO: array coulomb log
    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            fundamental_electron_collision_freq, insert_some_nans, insert_all_nans, {}
        )


class Test_fundamental_ion_collision_freq:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T_arr = np.array([1, 2]) * u.eV
        cls.n_arr = np.array([1e20, 2e20]) * u.cm**-3
        cls.ion = "p"
        cls.coulomb_log = 10

    # TODO: array coulomb log
    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            fundamental_ion_collision_freq, insert_some_nans, insert_all_nans, {}
        )
