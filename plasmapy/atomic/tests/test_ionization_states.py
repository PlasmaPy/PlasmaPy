import pytest
import numpy as np
import astropy.units as u
from ..ionization_states import IonizationState, IonizationStates
from ...utils import AtomicError, RunTestError, InvalidIsotopeError, run_test
from ...atomic import (
    atomic_number,
    mass_number,
    atomic_symbol,
    isotope_symbol,
    Particle,
)
import collections
import warnings


def has_attribute(attribute, tests_dict):
    cases = [test for test in tests_dict.keys() if attribute in tests_dict[test].keys()]
    if cases:
        return cases
    else:
        raise ValueError(f"No cases with attribute {attribute}")


tests = {

    'basic': {
        'inputs': {'H': [0.1, 0.9], 'helium': (0.2, 0.3, 0.5)},
        'T_e': 1e6 * u.K,
        'tol': 1e-15,
        'abundances': {'hydrogen': 1.0, 'He': 0.1},
    },

    'quantities': {
        'inputs': {'H': np.array([0.1, 0.9]) * u.m ** -3},
        'abundances': {'H': 1.0, 'He': 0.1}
    },

    'just H': {
        'inputs': {'H': [0.1, 0.9]},
    },

    'H acceptable error': {
        'inputs': {'H': [1.0, 1e-6]},
        'tol': 1e-5,
    },

    'n_H': {
        'inputs': {'H': [1, 0], 'He': [1, 0, 0]},
        'abundances': {'H': 1, 'He': 0.1},
        'n_H': 1e9 * u.cm **-3,
    },

    'T_e and n_H': {
        'inputs': {'H': [0.9, 0.1], 'He': [0.5, 0.3, 0.2]},
        'abundances': {'H': 1, 'He': 0.1},
        'T_e': 1e4 * u.K,
        'n_H': 1e15 * u.m ** -3,
    },

}

test_names2 = tests.keys()


class Test_IonizationStates:

    @classmethod
    def setup_class(cls):
        cls.instances = {}

    @pytest.mark.parametrize('test', test_names2)
    def test_instantiation(self, test):
        try:
            self.instances[test] = IonizationStates(**tests[test])
        except Exception as exc:
            raise AtomicError(
                f"Cannot create IonizationStates instance for "
                f"test='{test}'") from exc

    @pytest.mark.parametrize('test', test_names2)
    def test_keys(self, test):
        input_keys = tests[test]['inputs'].keys()
        particles = [Particle(input_key) for input_key in input_keys]
        expected_keys = [p.particle for p in particles]
        actual_keys = [key for key in self.instances[test].ionic_fractions.keys()]

        assert actual_keys == expected_keys, (
            f"For test='{test}', the following should be equal:\n"
            f"  actual_keys = {actual_keys}\n"
            f"expected_keys = {expected_keys}"
        )

    @pytest.mark.parametrize('test', has_attribute('abundances', tests))
    def test_abundances(self, test):
        try:
            actual_abundances = self.instances[test].abundances
        except Exception as exc:
            raise AttributeError("Unable to access abundances.") from exc

        elements = set(self.instances[test].elements)
        elements_from_abundances = set(actual_abundances.keys())

        if not elements.issubset(elements_from_abundances):
            raise RunTestError(
                f"The elements whose IonizationStates are being kept "
                f"track of ({elements}) are not a subset of the "
                f"elements whose abundances are being kept track of "
                f"({elements_from_abundances}) for test {test}."
            )

    @pytest.mark.parametrize('test', test_names2)
    def test_element_sorting(self, test):
        elements = self.instances[test].elements
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
            f"Elements/isotopes are not sorted for test='{test}':\n" 
            f"  before_sorting = {before_sorting}\n"
            f"   after_sorting = {after_sorting}\n"
            f"where above is (atomic_number, mass_number if isotope else 0)")

    @pytest.mark.parametrize('test', test_names2)
    def test_ionic_fractions(self, test):

        errmsg = ""

        elements = self.instances[test].elements
        input_keys = tests[test]["inputs"].keys()
        for element, input_key in zip(elements, input_keys):

            expected = np.array(tests[test]["inputs"][input_key])

            if isinstance(expected, u.Quantity):
                expected = np.array(expected.value / np.sum(expected.value))

            #if not isinstance(expected, np.ndarray)

            actual = self.instances[test].ionic_fractions[element]

            if not np.allclose(actual, expected):
                errmsg += (
                    f"\n\nThere is a discrepancy in ionic fractions for "
                    f"({test}, {element}, {input_key})\n"
                    f"  expected = {expected}\n"
                    f"    actual = {actual}")
            if not isinstance(actual, np.ndarray) or isinstance(actual, u.Quantity):
                raise AtomicError(
                    f"\n\nNot a numpy.ndarray: ({test}, {element})")

        if errmsg:
            raise AtomicError(errmsg)

#    @pytest.mark.parametrize('test', [test for test in tests if ])

    @pytest.mark.parametrize('test', test_names2)
    def test_getitem_element(self, test):
        """Test that __get_item__ returns an IonizationState instance"""
        instance = self.instances[test]

        for key in instance.elements:

            try:
                expected = instance.ionic_fractions[key]
            except Exception as exc:
                raise AtomicError(
                    f"Unable to get ionic_fractions for '{key}' in "
                    f"test='{test}'.") from exc

            try:
                actual = instance[key].ionic_fractions
            except Exception as exc:
                raise AtomicError(f"Unable to get item {key} in test={test}.")

            try:
                test_passed = np.allclose(expected, actual)
            except Exception:
                raise TypeError(
                    f"For test='{test}' and key='{key}', cannot "
                    f"compare expected ionic fractions of {expected} "
                    f"with the resulting ionic fractions of {actual}.") from None

            if not test_passed:
                raise AtomicError(
                    f"For test='{test}' and key='{key}', the expected "
                    f"ionic fractions of {expected} are not all equal "
                    f"to the resulting ionic fractions of {actual}.")

    @pytest.mark.parametrize('test', test_names2)
    def test_getitem_element_intcharge(self, test):
        instance = self.instances[test]
        for particle in instance.elements:
            for int_charge in range(0, atomic_number(particle) + 1):
                actual = instance[particle, int_charge].ionic_fraction
                expected = instance.ionic_fractions[particle][int_charge]
                # We only need to check if one is broken
                assert np.isclose(actual, expected), (
                    f"Indexing broken for:\n"
                    f"       test = '{test}'\n"
                    f"   particle = '{particle}'")

    @pytest.mark.parametrize('test', test_names2)
    def test_normalization(self, test):
        instance = self.instances[test]
        instance.normalize()
        not_normalized_elements = []
        for element in instance.elements:
            is_not_normalized = not np.isclose(
                np.sum(instance.ionic_fractions[element]),
                1,
                atol=1e-19,
                rtol=0,
            )

            if is_not_normalized:
                not_normalized_elements.append(element)

        if not_normalized_elements:
            raise AtomicError(
                f"In test='{test}', ionic fractions for the following "
                f"particles were not normalized: "
                f"{', '.join(not_normalized_elements)}."
            )



IE = collections.namedtuple("IE", ["inputs", "expected_exception"])

tests_for_exceptions = {
    'wrong type': IE({"inputs": None}, TypeError),
    'not normalized': IE({"inputs": {'He': [0.4, 0.5, 0.0]}, "tol": 1e-9}, AtomicError),
    'negative ionfrac': IE({"inputs": {'H': [-0.1, 1.1]}}, AtomicError),
    'ion': IE({"inputs": {'H': [0.1, 0.9], 'He+': [0.0, 0.9, 0.1]}}, AtomicError),
    'repeat elements': IE({"inputs": {'H': [0.1, 0.9], "hydrogen": [0.2, 0.8]}}, AtomicError),
    'isotope of element': IE({"inputs": {'H': [0.1, 0.9], "D": [0.2, 0.8]}}, AtomicError),

    'negative abundance': IE({
        "inputs": {"H": [0.1, 0.9], "He": [0.4, 0.5, 0.1]}, "abundances": {"H": 1, "He": -0.1},
    }, AtomicError),

    'imaginary abundance': IE({
        "inputs": {"H": [0.1, 0.9], "He": [0.4, 0.5, 0.1]}, "abundances": {"H": 1, "He": 0.1j},
    }, TypeError),

}


@pytest.mark.parametrize('test', tests_for_exceptions.keys())
def test_execeptions(test):
    """
    Test that appropriate exceptions are raised for inappropriate inputs
    to IonizationStates.
    """
    run_test(
        IonizationStates,
        kwargs=tests_for_exceptions[test].inputs,
        expected_outcome=tests_for_exceptions[test].expected_exception,
    )
