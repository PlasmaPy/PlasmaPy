import astropy.units as u
import itertools
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose
from numbers import Real
from typing import Dict

from plasmapy.particles import (
    atomic_number,
    IonicLevel,
    IonizationState,
    IonizationStateCollection,
    mass_number,
    Particle,
    particle_symbol,
    ParticleList,
)
from plasmapy.particles.exceptions import InvalidIsotopeError, ParticleError


def check_abundances_consistency(
    abundances: Dict[str, Real], log_abundances: Dict[str, Real]
):
    """
    Test that a set of abundances is consistent with a set of the base
    10 logarithm of abundances.
    """
    assert abundances.keys() == log_abundances.keys(), (
        f"Mismatch between keys from abundances and log_abundances.\n\n"
        f"abundances.keys():     {abundances.keys()}\n\n"
        f"log_abundances.keys(): {log_abundances.keys()}"
    )

    for element, abundance_from_abundances in abundances.items():
        abundance_from_log_abundances = 10 ** log_abundances[element]
        assert np.isclose(
            abundance_from_abundances, abundance_from_log_abundances
        ), "Mismatch between abundances and log_abundances."


def has_attribute(attribute, tests_dict):
    if cases := [test for test in tests_dict if attribute in tests_dict[test]]:
        return cases
    else:
        raise ValueError(f"No cases with attribute {attribute}")


tests = {
    "basic": {
        "inputs": {"H": [0.1, 0.9], "helium": (0.2, 0.3, 0.5)},
        "T_e": 1e6 * u.K,
        "tol": 1e-15,
        "abundances": {"hydrogen": 1.0, "He": 0.1},
        "kappa": 1.5001,
    },
    "quantities": {
        "inputs": {
            "H": np.array([10, 90]) * u.m**-3,
            "He": np.array([1, 9, 0]) * u.m**-3,
        }
    },
    "just H": {"inputs": {"H": [0.1, 0.9]}},
    "H acceptable error": {"inputs": {"H": [1.0, 1e-6]}, "tol": 1e-5},
    "n0": {
        "inputs": {"He": [1, 0, 0], "H": [1, 0]},
        "abundances": {"H": 1, "He": 0.1},
        "n0": 1e9 * u.cm**-3,
    },
    "T_e and n": {
        "inputs": {"H": [0.9, 0.1], "helium": [0.5, 0.3, 0.2]},
        "abundances": {"H": 1, "He": 0.1},
        "T_e": 1e4 * u.K,
        "n0": 1e15 * u.m**-3,
        "kappa": np.inf,
    },
    "log_abundances": {
        "inputs": {"H": [1, 0], "He": [1, 0, 0]},
        "log_abundances": {"H": 1, "He": 0},
        "n0": 1e9 * u.cm**-3,
    },
    "elements & isotopes": {
        "inputs": {
            "H": [0.9, 0.1],
            "He-3": [0.3, 0.7, 0.0],
            "He-4": [0.29, 0.69, 0.02],
        },
        "abundances": {"H": 1, "He-3": 1e-7, "He-4": 0.1},
        "n0": 1e12 * u.m**-3,
    },
    "ordered elements -> inputs": {"inputs": ["O", "C", "H", "Fe", "Ar"]},
    "mixed and unordered elements and isotopes": {
        "inputs": ("Og", "O", "H", "Fe-56", "He", "Li-7", "Li-6")
    },
    "number densities -> inputs": {
        "inputs": {
            "H": np.array([2, 3]) * u.m**-3,
            "He": np.array([5, 7, 11]) * u.m**-3,
        }
    },
    "number densities and n are both inputs": {
        "inputs": {"H": [0.1, 0.3] * u.cm**-3},
        "n0": 1e-5 * u.mm**-3,
    },
}

test_names = tests.keys()


@pytest.mark.slow
class TestIonizationStateCollection:
    @classmethod
    def setup_class(cls):
        cls.instances = {}

    @pytest.mark.parametrize("test_name", test_names)
    def test_instantiation(self, test_name):
        try:
            self.instances[test_name] = IonizationStateCollection(**tests[test_name])
        except Exception:
            pytest.fail(
                f"Cannot create IonizationStateCollection instance for test='{test_name}'"
            )

    @pytest.mark.parametrize("test_name", test_names)
    def test_no_exceptions_from_str(self, test_name):
        self.instances[test_name].__str__()

    @pytest.mark.parametrize("test_name", test_names)
    def test_no_exceptions_from_repr(self, test_name):
        self.instances[test_name].__repr__()

    @pytest.mark.parametrize("test_name", test_names)
    def test_no_exceptions_from_info(self, test_name):
        self.instances[test_name].summarize()

    @pytest.mark.parametrize("test_name", test_names)
    def test_simple_equality(self, test_name):
        """Test that __eq__ is not extremely broken."""
        a = IonizationStateCollection(**tests[test_name])
        b = IonizationStateCollection(**tests[test_name])
        assert a == a, "IonizationStateCollection instance does not equal itself."
        assert (
            a == b
        ), "IonizationStateCollection instance does not equal identical instance."

    @pytest.mark.parametrize(
        "test_name",
        [
            test_name
            for test_name in test_names
            if isinstance(tests[test_name]["inputs"], dict)
        ],
    )
    def test_that_particles_were_set_correctly(self, test_name):
        input_particles = tests[test_name]["inputs"].keys()
        particles = [Particle(input_particle) for input_particle in input_particles]
        expected_particles = {p.symbol for p in particles}
        actual_particles = set(self.instances[test_name].ionic_fractions)

        assert actual_particles == expected_particles, (
            f"For test='{test_name}', the following should be equal:\n"
            f"  actual_particles = {actual_particles}\n"
            f"expected_particles = {expected_particles}"
        )

    @pytest.mark.parametrize("test_name", has_attribute("abundances", tests))
    def test_that_abundances_kwarg_sets_abundances(self, test_name):
        try:
            actual_abundances = self.instances[test_name].abundances
        except Exception as exc:
            pytest.fail("Unable to access abundances.")

        elements = set(self.instances[test_name].base_particles)
        elements_from_abundances = set(actual_abundances.keys())

        if not elements.issubset(elements_from_abundances):
            pytest.fail(
                f"The elements whose IonizationStateCollection are being kept "
                f"track of ({elements}) are not a subset of the "
                f"elements whose abundances are being kept track of "
                f"({elements_from_abundances}) for test {test_name}."
            )

    @pytest.mark.parametrize("test_name", test_names)
    def test_that_elements_and_isotopes_are_sorted(self, test_name):
        elements = self.instances[test_name].base_particles
        before_sorting = []
        for element in elements:
            atomic_numb = atomic_number(element)
            try:
                mass_numb = mass_number(element)
            except InvalidIsotopeError:
                mass_numb = 0
            before_sorting.append((atomic_numb, mass_numb))
        after_sorting = sorted(before_sorting)

        assert before_sorting == after_sorting, (
            f"Elements/isotopes are not sorted for test='{test_name}':\n"
            f"  before_sorting = {before_sorting}\n"
            f"   after_sorting = {after_sorting}\n"
            f"where above is (atomic_number, mass_number if isotope else 0)"
        )

    @pytest.mark.parametrize("test_name", test_names)
    def test_that_ionic_fractions_are_set_correctly(self, test_name):

        errmsg = ""

        elements_actual = self.instances[test_name].base_particles
        inputs = tests[test_name]["inputs"]

        if isinstance(inputs, dict):
            input_keys = list(tests[test_name]["inputs"].keys())

            def sort_key(k):
                return atomic_number(k), mass_number(k) if Particle(k).isotope else 0

            input_keys = sorted(input_keys, key=sort_key)

            for element, input_key in zip(elements_actual, input_keys):
                expected = tests[test_name]["inputs"][input_key]

                if isinstance(expected, u.Quantity):
                    expected = np.array(expected.value / np.sum(expected.value))

                actual = self.instances[test_name].ionic_fractions[element]

                if not np.allclose(actual, expected):
                    errmsg += (
                        f"\n\nThere is a discrepancy in ionic fractions for "
                        f"({test_name}, {element}, {input_key})\n"
                        f"  expected = {expected}\n"
                        f"    actual = {actual}"
                    )

                if not isinstance(actual, np.ndarray) or isinstance(actual, u.Quantity):
                    raise ParticleError(
                        f"\n\nNot a numpy.ndarray: ({test_name}, {element})"
                    )
        else:
            elements_expected = {particle_symbol(element) for element in inputs}

            assert set(self.instances[test_name].base_particles) == elements_expected

            for element in elements_expected:
                assert all(np.isnan(self.instances[test_name].ionic_fractions[element]))
        if errmsg:
            pytest.fail(errmsg)

    @pytest.mark.parametrize("test_name", test_names)
    def test_getitem_element(self, test_name):
        """Test that __get_item__ returns an IonizationState instance"""
        instance = self.instances[test_name]

        for key in instance.base_particles:

            try:
                expected = instance.ionic_fractions[key]
            except Exception as exc:
                pytest.fail(
                    f"Unable to get ionic_fractions for '{key}' in test='{test_name}'."
                )

            try:
                actual = instance[key].ionic_fractions
            except Exception as exc:
                pytest(f"Unable to get item {key} in test={test_name}.")

            try:
                if all(np.isnan(expected)):
                    test_passed = True
                else:
                    test_passed = np.allclose(expected, actual)
            except Exception:
                raise TypeError(
                    f"For test='{test_name}' and key='{key}', cannot "
                    f"compare expected ionic fractions of {expected} "
                    f"with the resulting ionic fractions of {actual}."
                ) from None

            if not test_passed:
                pytest.fail(
                    f"For test='{test_name}' and key='{key}', the expected "
                    f"ionic fractions of {expected} are not all equal "
                    f"to the resulting ionic fractions of {actual}."
                )

    @pytest.mark.parametrize("test_name", test_names)
    def test_getitem_element_intcharge(self, test_name):
        instance = self.instances[test_name]
        for particle in instance.base_particles:
            for int_charge in range(atomic_number(particle) + 1):
                actual = instance[particle, int_charge].ionic_fraction
                expected = instance.ionic_fractions[particle][int_charge]
                # We only need to check if one is broken
            if not np.isnan(actual) and np.isnan(expected):
                assert np.isclose(actual, expected), (
                    f"Indexing broken for:\n"
                    f"       test = '{test_name}'\n"
                    f"   particle = '{particle}'"
                )

    @pytest.mark.parametrize(
        "test_name",
        [
            test_name
            for test_name in test_names
            if isinstance(tests[test_name]["inputs"], dict)
        ],
    )
    def test_normalization(self, test_name):
        instance = self.instances[test_name]
        instance.normalize()
        not_normalized_elements = []
        for element in instance.base_particles:
            is_not_normalized = not np.isclose(
                np.sum(instance.ionic_fractions[element]), 1, atol=1e-19, rtol=0
            )

            if is_not_normalized:
                not_normalized_elements.append(element)

        if not_normalized_elements:
            pytest.fail(
                f"In test = '{test_name}', ionic fractions for the "
                f"following particles were not normalized: "
                f"{', '.join(not_normalized_elements)}."
            )


def test_abundances_consistency():
    """Test that ``abundances`` and ``log_abundances`` are consistent."""

    inputs = {"H": [1, 0], "He": [1, 0, 0]}
    abundances = {"H": 1.0, "He": 0.1}
    elements = abundances.keys()

    log_abundances = {element: np.log10(abundances[element]) for element in elements}

    instance_nolog = IonizationStateCollection(inputs, abundances=abundances)
    instance_log = IonizationStateCollection(inputs, log_abundances=log_abundances)

    for element in elements:
        assert np.allclose(
            instance_log.abundances[element], instance_nolog.abundances[element]
        ), "abundances not consistent."

    for element in elements:
        assert np.allclose(
            instance_log.log_abundances[element], instance_nolog.log_abundances[element]
        ), "log_abundances not consistent."


class TestIonizationStateCollectionItemAssignment:
    """
    Test IonizationStateCollection.__setitem__ and exceptions.
    """

    @classmethod
    def setup_class(cls):
        cls.states = IonizationStateCollection(
            {"H": [0.9, 0.1], "He": [0.5, 0.4999, 1e-4]}
        )

    @pytest.mark.parametrize(
        "element, new_states",
        [
            ("H", [np.nan, np.nan]),
            ("He", [np.nan, np.nan, np.nan]),
            ("H", [0.1, 0.9]),
            ("He", [0.89, 0.1, 0.01]),
        ],
    )
    def test_setitem(self, element, new_states):
        """Test item assignment in an IonizationStateCollection instance."""
        try:
            self.states[element] = new_states
        except Exception:
            pytest.fail(
                "Unable to change ionic fractions for an IonizationStateCollection instance."
            )
        resulting_states = self.states[element].ionic_fractions

        assert np.any(
            [
                np.allclose(resulting_states, new_states),
                np.all(np.isnan(resulting_states)) and np.all(np.isnan(new_states)),
            ]
        )

    @pytest.mark.parametrize(
        "base_particle, new_states, expected_exception",
        [
            ("H", (0, 0.9), ValueError),
            ("H", (-0.1, 1.1), ValueError),
            ("H", (0.0, 1.0, 0.0), ValueError),
            ("Li", (0.0, 1.0, 0.0, 0.0), KeyError),
            ("sdfasd", (0, 1), KeyError),
            (KeyError, KeyError, KeyError),
        ],
    )
    def test_setitem_errors(self, base_particle, new_states, expected_exception):
        with pytest.raises(expected_exception):
            self.states[base_particle] = new_states


class TestIonizationStateCollectionDensities:
    @classmethod
    def setup_class(cls):

        cls.initial_ionfracs = {
            "H": np.array([0.87, 0.13]),
            "He": np.array([0.24, 0.37, 0.39]),
        }
        cls.abundances = {"H": 1.0, "He": 0.0835}
        cls.n = 10 * u.m**-3

        cls.expected_densities = {
            "H": np.array([8.7, 1.3]) * u.m**-3,
            "He": np.array([0.2004, 0.30895, 0.32565]) * u.m**-3,
        }

        cls.expected_electron_density = 2.26025 * u.m**-3
        cls.states = IonizationStateCollection(
            cls.initial_ionfracs, abundances=cls.abundances, n0=cls.n
        )

    def test_electron_density(self):
        assert np.isclose(
            self.states.n_e.value, self.expected_electron_density.value
        ), (
            "Mismatch in electron density calculation:\n"
            f"Calculated = {self.states.n_e}\n"
            f"Expected   = {self.expected_electron_density}"
        )

    @pytest.mark.parametrize("elem", ["H", "He"])
    def test_number_densities(self, elem):
        assert np.allclose(
            self.states.number_densities[elem].value,
            self.expected_densities[elem].value,
        ), (
            f"Mismatch in number densities for {elem}\n"
            f"Calculated = {self.states.number_densities[elem]}\n"
            f"Expected   = {self.expected_electron_density}"
        )


class TestIonizationStateCollectionAttributes:
    @classmethod
    def setup_class(cls):
        cls.elements = ["H", "He", "Li", "Fe"]
        cls.instance = IonizationStateCollection(cls.elements)
        cls.new_n = 5.153 * u.cm**-3

    @pytest.mark.parametrize("uninitialized_attribute", ["T_e", "n0", "n_e"])
    def test_attribute_defaults_to_nan(self, uninitialized_attribute):
        command = f"self.instance.{uninitialized_attribute}"
        default_value = eval(command)
        assert np.isnan(default_value), (
            f"{uninitialized_attribute} does not default to nan but "
            f"instead defaults to {default_value}."
        )

    def test_kappa_defaults_to_inf(self):
        assert np.isinf(
            self.instance.kappa
        ), "kappa does not default to a value of inf."

    @pytest.mark.parametrize(
        "uninitialized_attribute", ["number_densities", "ionic_fractions"]
    )
    def test_attribute_defaults_to_dict_of_nans(self, uninitialized_attribute):
        command = f"self.instance.{uninitialized_attribute}"
        default_value = eval(command)
        assert (
            list(default_value.keys()) == self.elements
        ), "Incorrect base particle keys."
        for element in self.elements:
            assert (
                len(default_value[element]) == atomic_number(element) + 1
            ), f"Incorrect number of ionization levels for {element}."
            assert np.all(np.isnan(default_value[element])), (
                f"The values do not default to an array of nans for " f"{element}."
            )

    @pytest.mark.parametrize(
        "uninitialized_attribute", ["abundances", "log_abundances"]
    )
    def test_abundances_default_to_nans(self, uninitialized_attribute):
        command = f"self.instance.{uninitialized_attribute}"
        default_value = eval(command)
        for element in self.elements:
            assert isinstance(default_value[element], Real)
            assert np.isnan(default_value[element])

    @pytest.mark.parametrize(
        "attribute, invalid_value, expected_exception",
        [
            ("T_e", "5 * u.m", u.UnitsError),
            ("T_e", "-1 * u.K", ParticleError),
            ("n0", "5 * u.m", u.UnitsError),
            ("n0", "-1 * u.m ** -3", ParticleError),
            (
                "ionic_fractions",
                {"H": [0.3, 0.7], "He": [-0.1, 0.4, 0.7]},
                ParticleError,
            ),
            (
                "ionic_fractions",
                {"H": [0.3, 0.7], "He": [1.01, 0.0, 0.7]},
                ParticleError,
            ),
            (
                "ionic_fractions",
                {"H": [0.3, 0.6], "He": [1.0, 0.0, 0.0]},
                ParticleError,
            ),
            ("ionic_fractions", {"H": [1.0, 0.0]}, ParticleError),
        ],
    )
    def test_attribute_exceptions(self, attribute, invalid_value, expected_exception):

        command = f"self.instance.{attribute} = {invalid_value}"
        errmsg = f"No {expected_exception} was raised for command\n\n: {command}"

        with pytest.raises(expected_exception):
            exec(command)
            pytest.fail(errmsg)

    def test_setting_ionic_fractions_for_single_element(self):
        """
        Test that __setitem__ correctly sets new ionic fractions when
        used for just H, while not changing the ionic fractions for He
        from the uninitialized default of an array of nans of length 3.
        """
        self.new_fractions = [0.3, 0.7]
        self.instance["H"] = self.new_fractions
        resulting_fractions = self.instance.ionic_fractions["H"]
        assert np.allclose(
            self.new_fractions, resulting_fractions
        ), "Ionic fractions for H not set using __setitem__."
        assert "He" in self.instance.ionic_fractions, (
            "He is missing in ionic_fractions after __setitem__ was "
            "used to set H ionic fractions."
        )
        assert np.all(np.isnan(self.instance.ionic_fractions["He"])), (
            "He ionic fractions are not all nans after __setitem__ "
            "was used to set H ionic fractions."
        )

    @pytest.mark.parametrize(
        "key, invalid_fracs, expected_exception",
        [
            ("H", [-0.01, 1.01], ValueError),
            ("H", [0.4, 0.5], ValueError),
            ("H", [0.5, 0.5, 0.0], ValueError),
            ("He", [0.5, 0.5], ValueError),
            ("He", [0.1, 0.9, 0.0, 0.0], ValueError),
            ("He", [0.9, 0.1, 0.0] * u.m**-2, ValueError),
            ("He", [-0.01, 0.99, 0.02], ValueError),
            ("He", [1.01, -0.02, 0.01], ValueError),
        ],
    )
    def test_setting_invalid_ionfracs(self, key, invalid_fracs, expected_exception):
        errmsg = (
            f"No {expected_exception} is raised when trying to assign "
            f"{invalid_fracs} to {key} in an IonizationStateCollection instance."
        )
        with pytest.raises(expected_exception):
            self.instance[key] = invalid_fracs
            pytest.fail(errmsg)

    def test_setting_incomplete_abundances(self):
        new_abundances = {"H": 1, "He": 0.1, "Fe": 1e-5, "Au": 1e-8}  # missing lithium
        with pytest.raises(ParticleError):
            self.instance.abundances = new_abundances

    def test_setting_abundances(self):
        new_abundances = {"H": 1, "He": 0.1, "Li": 1e-4, "Fe": 1e-5, "Au": 1e-8}

        log_new_abundances = {
            element: np.log10(new_abundances[element]) for element in new_abundances
        }

        try:
            self.instance.abundances = new_abundances
        except Exception:
            pytest(f"Could not set abundances to {new_abundances}.")
        else:
            check_abundances_consistency(
                self.instance.abundances, self.instance.log_abundances
            )

        try:
            self.instance.log_abundances = log_new_abundances
        except Exception:
            pytest.fail(f"Could not set log_abundances to {log_new_abundances}.")
        else:
            check_abundances_consistency(
                self.instance.abundances, self.instance.log_abundances
            )

    @pytest.mark.parametrize(
        "invalid_indices",
        [
            (1, 2, 3),
            "C",
            "H-1",
            ("Fe", -1),
            ("Fe", 27),
            ("He", -1),
            ("He", 3),
            ("Fe", slice(3, 7)),
        ],
    )
    def test_invalid_indices(self, invalid_indices):
        with pytest.raises(IndexError):
            self.instance[invalid_indices]

    @pytest.mark.parametrize("index", ["H", "Fe"])
    def test_getitem_one_index(self, index):
        instance = self.instance
        result = instance[index]

        if np.all(np.isnan(instance.number_densities[index])):
            inputs = instance.ionic_fractions[index]
        else:
            inputs = instance.number_densities[index]

        expected = IonizationState(
            index, inputs, T_e=instance.T_e, kappa=instance.kappa
        )

        assert isinstance(result, IonizationState)
        assert result == expected

    @pytest.mark.parametrize("indices", [("H", 1), ("Fe", 6)])
    def test_getitem_two_indices(self, indices):
        instance = self.instance
        result = instance[indices]

        particle = indices[0]
        charge_number = indices[1]

        assert isinstance(result, IonicLevel)
        assert result.charge_number == charge_number

        expected_ionic_fraction = instance.ionic_fractions[particle][charge_number]

        assert np.any(
            [
                np.isclose(result.ionic_fraction, expected_ionic_fraction),
                np.isnan(result.ionic_fraction) and np.isnan(expected_ionic_fraction),
            ]
        )

        assert result.ionic_symbol == particle_symbol(particle, Z=charge_number)

    def test_setting_n(self):
        try:
            self.instance.n0 = self.new_n
        except Exception:
            pytest.fail("Unable to set number density scaling factor attribute")
        if not u.quantity.allclose(self.instance.n0, self.new_n):
            pytest.fail("Number density scaling factor was not set correctly.")
        if self.instance.n0.unit != u.m**-3:
            pytest.fail("Incorrect units for new number density.")

    def test_resetting_valid_densities(self):
        """
        Test that item assignment can be used to set number densities
        that preserve the total element number density.
        """

        element = "H"
        valid_ionic_fractions = [0.54, 0.46]
        original_n_elem = np.sum(self.instance.number_densities[element])
        valid_number_densities = valid_ionic_fractions * original_n_elem

        try:
            self.instance[element] = valid_number_densities
        except Exception:
            pytest.fail("Unable to set valid number densities using item assignment.")

        assert u.quantity.allclose(
            self.instance.ionic_fractions[element], valid_ionic_fractions
        ), "Item assignment of valid number densities did not yield correct ionic fractions."

        assert u.quantity.allclose(
            self.instance.number_densities[element], valid_number_densities
        ), "Item assignment of valid number densities did not yield correct number densities."

    def test_resetting_invalid_densities(self):
        """
        Test that item assignment with number densities that would
        change the total element number density raises an exception.
        """
        element = "H"
        original_n_elem = np.sum(self.instance.number_densities[element])
        invalid_number_densities = np.array([1.0001, 0]) * original_n_elem
        with pytest.raises(ValueError):
            self.instance[element] = invalid_number_densities

    def test_elemental_abundances_not_quantities(self):
        for element in self.instance.base_particles:
            assert not isinstance(self.instance.abundances[element], u.Quantity)

    @pytest.mark.parametrize("element", ["H", "He", "Fe"])
    def test_ionic_fractions_not_quantities(self, element):
        ionic_fractions = self.instance.ionic_fractions[element]
        if isinstance(ionic_fractions, u.Quantity):
            pytest.fail(
                f"The ionic fractions of {element} are a Quantity but should not be."
            )

    def test_that_iron_ionic_fractions_are_still_undefined(self):
        assert "Fe" in self.instance.ionic_fractions
        iron_fractions = self.instance.ionic_fractions["Fe"]
        assert len(iron_fractions) == atomic_number("Fe") + 1
        assert np.all(np.isnan(iron_fractions))

    def test_base_particles(self):
        """
        Test that the original base particles remain as the base
        particles after performing a bunch of operations that should not
        change them.
        """
        assert self.instance.base_particles == self.elements

    def test_base_particles_equal_ionic_fraction_particles(self):
        assert self.instance.base_particles == list(
            self.instance.ionic_fractions.keys()
        )


class TestIonizationStateCollectionDensityEqualities:
    """
    Test that IonizationStateCollection instances are equal or not equal to each
    other as they should be for different combinations of inputs
    related to ionic_fractions, number densities, and abundances.
    """

    @classmethod
    def setup_class(cls):

        # Create arguments to IonizationStateCollection that are all consistent
        # with each other.

        cls.ionic_fractions = {"H": [0.9, 0.1], "He": [0.99, 0.01, 0.0]}
        cls.abundances = {"H": 1, "He": 0.08}
        cls.n = 5.3 * u.m**-3
        cls.number_densities = {
            element: cls.ionic_fractions[element] * cls.n * cls.abundances[element]
            for element in cls.ionic_fractions
        }

        # The keys that begin with 'ndens' have enough information to
        # yield the number_densities attribute, whereas the keys that
        # begin with "no_ndens" do not.

        cls.dict_of_kwargs = {
            "ndens1": {
                "inputs": cls.ionic_fractions,
                "abundances": cls.abundances,
                "n0": cls.n,
            },
            "ndens2": {"inputs": cls.number_densities},
            "no_ndens3": {"inputs": cls.ionic_fractions},
            "no_ndens4": {"inputs": cls.ionic_fractions, "abundances": cls.abundances},
            "no_ndens5": {"inputs": cls.ionic_fractions, "n0": cls.n},
        }

        cls.instances = {
            key: IonizationStateCollection(**cls.dict_of_kwargs[key])
            for key in cls.dict_of_kwargs
        }

    @pytest.mark.parametrize("test_key", ["ndens1", "ndens2"])
    def test_number_densities_defined(self, test_key):
        number_densities = self.instances[test_key].number_densities
        for base_particle in self.instances[test_key].base_particles:
            assert not np.any(np.isnan(number_densities[base_particle])), (
                f"Test {test_key} should have number densities "
                f"defined, but doesn't."
            )

    @pytest.mark.parametrize("test_key", ["no_ndens3", "no_ndens4", "no_ndens5"])
    def test_number_densities_undefined(self, test_key):
        number_densities = self.instances[test_key].number_densities
        for base_particle in self.instances[test_key].base_particles:
            assert np.all(np.isnan(number_densities[base_particle])), (
                f"Test {test_key} should not have number densities "
                f"defined, but does."
            )

    @pytest.mark.parametrize(
        "this, that",
        itertools.product(
            ["ndens1", "ndens2", "no_ndens3", "no_ndens4", "no_ndens5"], repeat=2
        ),
    )
    def test_equality(self, this, that):
        """
        Test that the IonizationStateCollection instances that should provide
        ``number_densities`` are all equal to each other.  Test that the
        instances that should not provide ``number_densities`` are all
        equal to each other.  Test that each instance that should
        provide ``number_densities`` is not equal to each instance that
        should not provide ``number_densities``.
        """
        expect_equality = this[:4] == that[:4]
        are_equal = self.instances[this] == self.instances[that]
        if expect_equality != are_equal:
            print(f"{this} kwargs:\n {self.dict_of_kwargs[this]}\n")
            self.instances[this].summarize()
            print()
            print(f"{that} kwargs:\n {self.dict_of_kwargs[that]}\n")
            self.instances[that].summarize()
            descriptor = "equal" if expect_equality else "unequal"
            pytest.fail(
                f"Cases {this} and {that} should be {descriptor} but " f"are not."
            )


def test_number_density_assignment():
    instance = IonizationStateCollection(["H", "He"])
    number_densities = [2, 3, 5] * u.m**-3
    instance["He"] = number_densities


def test_len():
    ionization_states = IonizationStateCollection(["H", "He"])
    assert len(ionization_states) == 2


def test_iteration_with_nested_iterator():
    ionization_states = IonizationStateCollection(["H", "He"])

    i = 0
    for ionization_state1 in ionization_states:
        assert isinstance(ionization_state1, IonizationState)
        for ionization_state2 in ionization_states:
            assert isinstance(ionization_state2, IonizationState)
            i += 1
    assert i == 4


@pytest.mark.xfail()
def test_two_isotopes_of_same_element():
    instance = IonizationStateCollection(["H-1", "D"])


example_ionic_fractions = [
    ("H-1", [0, 1]),
    ("H-1", [0.2, 0.8]),
    ("He-4", [0.2, 0.3, 0.5]),
]


physical_properties = ["charge", "mass"]


@pytest.mark.parametrize("include_neutrals", [True, False])
@pytest.mark.parametrize("use_rms_mass", [False, True])
@pytest.mark.parametrize("use_rms_charge", [False, True])
@pytest.mark.parametrize("base_particle, ionic_fractions", example_ionic_fractions)
@pytest.mark.parametrize("physical_type", physical_properties)
def test_average_ion_consistency(
    base_particle,
    ionic_fractions,
    include_neutrals,
    use_rms_mass,
    use_rms_charge,
    physical_type,
):
    """
    Make sure that the average ions returned from equivalent `IonizationState`
    and `IonizationStateCollection` instances are consistent with each other.
    """
    ionization_state = IonizationState(base_particle, ionic_fractions)
    ionization_state_collection = IonizationStateCollection(
        {base_particle: ionic_fractions}
    )

    options = {
        "include_neutrals": include_neutrals,
        "use_rms_charge": use_rms_charge,
        "use_rms_mass": use_rms_mass,
    }

    average_ion_from_ionization_state = ionization_state.average_ion(**options)
    average_ion_from_ionization_state_collection = (
        ionization_state_collection.average_ion(**options)
    )

    quantity_from_ion_state = getattr(average_ion_from_ionization_state, physical_type)
    quantity_from_ion_collection = getattr(
        average_ion_from_ionization_state_collection, physical_type
    )

    assert_quantity_allclose(
        quantity_from_ion_state, quantity_from_ion_collection, rtol=1e-10
    )


@pytest.mark.parametrize("physical_property", physical_properties)
@pytest.mark.parametrize("use_rms", [True, False])
@pytest.mark.parametrize("include_neutrals", [True, False])
def test_comparison_to_equivalent_particle_list(
    physical_property, use_rms, include_neutrals
):
    """
    Test that `IonizationState.average_ion` gives consistent results with
    `ParticleList.average_particle` when the ratios of different particles
    is the same between the `IonizationState` and the `ParticleList`.
    """

    neutrals = 3 * ["H-1 0+"] + 2 * ["He-4 0+"] if include_neutrals else []
    ions = 2 * ["p+"] + 3 * ["He-4 1+"] + 5 * ["α"]
    particles = ParticleList(neutrals + ions)

    ionic_fractions = {
        "H-1": [0.6, 0.4],
        "He-4": [0.2, 0.3, 0.5],
    }

    abundances = {"H-1": 1, "He-4": 2}
    ionization_state_collection = IonizationStateCollection(
        ionic_fractions, abundances=abundances
    )

    kwarg = {f"use_rms_{physical_property}": True}
    expected_average_particle = particles.average_particle(**kwarg)
    expected_average_quantity = getattr(expected_average_particle, physical_property)

    actual_average_particle = ionization_state_collection.average_ion(
        include_neutrals=include_neutrals, **kwarg
    )
    actual_average_quantity = getattr(actual_average_particle, physical_property)

    assert_quantity_allclose(actual_average_quantity, expected_average_quantity)


def test_average_particle_exception():
    """
    Test that `IonizationStateCollection.average_ion` raises the
    appropriate exception when abundances are undefined and there is more
    than one base element or isotope.
    """
    ionization_states = IonizationStateCollection({"H": [1, 0], "He": [1, 0, 0]})

    with pytest.raises(ParticleError):
        ionization_states.average_ion()


def test_inequality_with_different_type():
    ionization_states = IonizationStateCollection({"H": [1, 0], "He": [1, 0, 0]})
    assert ionization_states != "different type"


def test_inequality_with_different_base_particles():
    instance1 = IonizationStateCollection({"H": [1, 0], "He": [1, 0, 0]})
    instance2 = IonizationStateCollection({"H": [1, 0], "Li": [1, 0, 0, 0]})
    assert instance1 != instance2
