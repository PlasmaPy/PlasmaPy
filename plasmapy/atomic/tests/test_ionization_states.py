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
    particle_symbol,
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

    'n': {
        'inputs': {'H': [1, 0], 'He': [1, 0, 0]},
        'abundances': {'H': 1, 'He': 0.1},
        'n': 1e9 * u.cm **-3,
    },

    'T_e and n': {
        'inputs': {'H': [0.9, 0.1], 'He': [0.5, 0.3, 0.2]},
        'abundances': {'H': 1, 'He': 0.1},
        'T_e': 1e4 * u.K,
        'n': 1e15 * u.m ** -3,
    },

    'log_abundances': {
        'inputs': {'H': [1, 0], 'He': [1, 0, 0]},
        'log_abundances': {'H': 1, 'He': 0},
        'n': 1e9 * u.cm ** -3,
    },

    'elements & isotopes': {
        'inputs': {'H': [0.9, 0.1], 'He-3': [0.3, 0.7, 0.0], 'He-4': [0.29, 0.69, 0.02]},
        'abundances': {'H': 1, 'He-3': 1e-7, 'He-4': 0.1},
        'n': 1e12 * u.m ** -3,
    },

    'just elements': {
        'inputs': ['H', 'He'],
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

    @pytest.mark.parametrize(
        'test',
        [name for name in test_names2 if isinstance(tests[name]['inputs'], dict)],
    )
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

        elements_actual = self.instances[test].elements
        inputs = tests[test]["inputs"]

        if isinstance(inputs, dict):

            input_keys = tests[test]["inputs"].keys()
            for element, input_key in zip(elements_actual, input_keys):

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

        else:
            elements_expected = {particle_symbol(element) for element in inputs}

            assert set(self.instances[test].elements) == elements_expected

            for element in elements_expected:
                assert all(np.isnan(self.instances[test].ionic_fractions[element]))

        if errmsg:
            raise AtomicError(errmsg)

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
                if all(np.isnan(expected)):
                    test_passed=True
                else:
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
#                if not (all(np.isnan(actual)) and all(np.isnan(expected))):
            if np.isnan(actual) and np.isnan(expected):
                continue
            else:
                assert np.isclose(actual, expected), (
                    f"Indexing broken for:\n"
                    f"       test = '{test}'\n"
                    f"   particle = '{particle}'")

    @pytest.mark.parametrize(
        'test',
        [name for name in test_names2 if isinstance(tests[name]['inputs'], dict)]
    )
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


def test_IonizationStates_abundances():
    """Test that abundances and log_abundances are consistent."""

    inputs = {'H': [1, 0], 'He': [1, 0, 0]}
    abundances = {'H': 1.0, 'He': 0.1}
    elements = abundances.keys()


    log_abundances = {element: np.log10(abundances[element]) for element in elements}

    instance_nolog = IonizationStates(inputs, abundances=abundances)
    instance_log = IonizationStates(inputs, log_abundances=log_abundances)

    for element in elements:
        assert np.allclose(
            instance_log.abundances[element],
            instance_nolog.abundances[element],
        ), 'abundances not consistent.'

    for element in elements:
        assert np.allclose(
            instance_log.log_abundances[element],
            instance_nolog.log_abundances[element],
        ), 'log_abundances not consistent.'


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


def test_setitem():
    states = IonizationStates({'H': [0.9, 0.1], 'He': [0.5, 0.4999, 1e-4]})

    new_states = [0.0, 1.0]
    states['H'] = new_states
    assert np.allclose(states['H'].ionic_fractions, new_states)

@pytest.mark.parametrize(
    'new_states,expected_exception',
    [
        ((0, 0.9), AtomicError),
        ((-0.1, 1.1), AtomicError),
        ((0.0, 1.0, 0.0), AtomicError),
     ]
)
def test_setitem_errors(new_states, expected_exception):
    states = IonizationStates({'H': [0.9, 0.1], 'He': [0.5, 0.4999, 1e-4]})
    with pytest.raises(expected_exception):
        states['H'] = new_states


class Test_IonizationStates:

    @classmethod
    def setup_class(cls):

        cls.initial_ionfracs = {'H': np.array([0.87, 0.13]), 'He': np.array([0.24, 0.37, 0.39])}
        cls.abundances = {'H': 1.0, 'He': 0.0835}
        cls.n = 10 * u.m ** -3

        cls.expected_densities = {
            'H': np.array([8.7, 1.3]) * u.m ** -3,
            'He': np.array([0.2004 , 0.30895, 0.32565]) * u.m ** -3
        }

        cls.expected_electron_density = 2.26025 * u.m ** -3
        cls.states = IonizationStates(cls.initial_ionfracs, abundances=cls.abundances, n=cls.n)

    def test_electron_density(self):
        assert np.isclose(self.states.n_e.value, self.expected_electron_density.value), (
            'Mismatch in electron density calculation:\n'
            f'Calculated = {self.states.n_e}\n'
            f'Expected   = {self.expected_electron_density}')

    @pytest.mark.parametrize('elem', ['H', 'He'])
    def test_number_densities(self, elem):
        assert np.allclose(
            self.states.number_densities[elem].value,
            self.expected_densities[elem].value), (
            f"Mismatch in number densities for {elem}\n"
            f"Calculated = {self.states.number_densities[elem]}\n"
            f"Expected   = {self.expected_electron_density}")



