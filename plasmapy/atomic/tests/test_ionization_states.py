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

test_cases = {
    'Li': {
        'particle': 'Li',
        'ionic_fractions': np.array([0.4, 0.3, 0.2, 0.1]),
        'tol': 1e-15,
    },

    'Li ground state': {
        'particle': 'Li',
        'ionic_fractions': np.array([1, 0, 0, 0], dtype=np.int64),
        'tol': 1e-15,
    },

    'H': {
        'particle': 'H',
        'ionic_fractions': [0.6, 0.4],
        'tol': 1e-8,
    },

    'H acceptable error': {
        'particle': 'H',
        'ionic_fractions': [0.6, 0.400_000_001],
        'tol': 1e-8,
    },

    'D': {
        'particle': 'deuterium',
        'ionic_fractions': [0.7, 0.3],
        'tol': 1e-15,
    },

    'He': {
        'particle': 'He',
        'ionic_fractions': [0.5, 0.3, 0.2],
        'n_elem': 1e20 * u.m ** -3,
    },

}

test_names = test_cases.keys()


class Test_IonizationState:
    """Test instances of IonizationState."""

    @classmethod
    def setup_class(cls):
        "Set up the class and test instantantiation."
        cls.instances = {}
        for test_name in test_names:
            try:
                cls.instances[test_name] = IonizationState(**test_cases[test_name])
            except Exception as exc:
                raise RunTestError(
                    f"Unable to create IonizationState instance for "
                    f"test case {test_name}.")

    @pytest.mark.parametrize('test_name', test_names)
    def test_ionic_fractions(self, test_name):
        instance = self.instances[test_name]
        assert np.allclose(instance.ionic_fractions, test_cases[test_name]['ionic_fractions'])

    def test_equality1(self):
        assert self.instances['Li'] == self.instances['Li'], \
            "Identical IonizationState instances are not equal."

    def test_equality2(self):
        assert self.instances['H'] == self.instances['H acceptable error'], \
            ("Two IonizationState instances that are approximately the "
             "same to within the tolerance are not testing as equal.")

    def test_inequality(self):
        assert self.instances['Li ground state'] != self.instances['Li'], \
            "Different IonizationState instances are equal."

    def test_equality_error(self):
        """
        Test that comparisons of IonizationState instances for
        different elements
        """
        with pytest.raises(AtomicError):
            self.instances['Li'] == self.instances['H']

    @pytest.mark.parametrize('test_name', test_names)
    def test_iteration(self, test_name: str):
        """Test that IonizationState instances iterate impeccably."""
        try:
            states = [state for state in self.instances[test_name]]
        except Exception:
            raise AtomicError(
                f"Unable to perform iteration for {test_name}.")

        try:
            integer_charges = [state.integer_charge for state in states]
            ionic_fractions = np.array([state.ionic_fraction for state in states])
            ionic_symbols = [state.ionic_symbol for state in states]
        except Exception:
            raise AtomicError("An attribute may be misnamed or missing.")

        try:
            base_symbol = isotope_symbol(ionic_symbols[0])
        except InvalidIsotopeError:
            base_symbol = atomic_symbol(ionic_symbols[0])
        finally:
            atomic_numb = atomic_number(ionic_symbols[1])

        errors = []

        expected_charges = np.arange(atomic_numb + 1)
        if not np.all(integer_charges == expected_charges):
            errors.append(
                f"The resulting integer charges are {integer_charges}, "
                f"which are not equal to the expected integer charges, "
                f"which are {expected_charges}.")

        expected_fracs = test_cases[test_name]['ionic_fractions']
        if not np.allclose(ionic_fractions, expected_fracs):
            errors.append(
                f"The resulting ionic fractions are {ionic_fractions}, "
                f"which are not equal to the expected ionic fractions "
                f"of {expected_fracs}.")

        expected_particles = [Particle(base_symbol, Z=charge) for charge in integer_charges]
        expected_symbols = [particle.ionic_symbol for particle in expected_particles]
        if not ionic_symbols == expected_symbols:
            errors.append(
                f"The resulting ionic symbols are {ionic_symbols}, "
                f"which are not equal to the expected ionic symbols of "
                f"{expected_symbols}.")

        if errors:
            errors.insert(0, (
                f"The test of IonizationState named '{test_name}' has "
                f"resulted in the following errors when attempting to "
                f"iterate."))
            errmsg = " ".join(errors)
            raise AtomicError(errmsg)

    def test_slicing1(self):
        assert np.allclose(
            self.instances['Li'][1:3].ionic_fraction,
            test_cases['Li']['ionic_fractions'][1:3]
        )

    def test_slicing2(self):
        assert np.allclose(
            self.instances['Li'][1:4:2].ionic_fraction,
            test_cases['Li']['ionic_fractions'][1:4:2]
        )

    @pytest.mark.parametrize('index', [-1, 4, 'Li'])
    def test_indexing_error(self, index):
        with pytest.raises(AtomicError):
            self.instances['Li'][index]

    def test_normalization(self):
        H = self.instances['H acceptable error']
        assert not H.is_normalized(tol=1e-15)
        H.normalize()
        assert H.is_normalized(tol=1e-15)

    @pytest.mark.parametrize('test_name', test_names)
    def test_identifications(self, test_name):
        """
        Test that the identification attributes for test IonizationState
        instances match the expected values from the Particle instance.
        """

        Identifications = collections.namedtuple(
            "Identifications",
            ["element", "isotope", "base_particle", "atomic_number"],
        )

        expected_identifications = Identifications(
            self.instances[test_name].element,
            self.instances[test_name].isotope,
            self.instances[test_name].base_particle,
            self.instances[test_name].atomic_number,
        )

        expected_element = self.instances[test_name]._particle.element
        expected_isotope = self.instances[test_name]._particle.isotope
        expected_atomic_number = self.instances[test_name]._particle.atomic_number

        resulting_identifications = Identifications(
            expected_element,
            expected_isotope,
            expected_isotope if expected_isotope else expected_element,
            expected_atomic_number,
        )

        assert resulting_identifications == expected_identifications, (
            f"For IonizationState test {test_name}, the resulting "
            f"identifications of {resulting_identifications} differ "
            f"from the expected identifications of "
            f"{expected_identifications}."
        )

    @pytest.mark.parametrize('tol', [-1e-16, 1.0000001])
    def test_invalid_tolerances(self, tol):
        """Test that invalid tolerances raise appropriate errors."""
        test_name = "Li"
        instance = self.instances[test_name]
        with pytest.raises(ValueError):
            instance.tol = tol

    @pytest.mark.parametrize('test_name', test_cases.keys())
    def test_particles(self, test_name):
        """
        Test that IonizationState returns the correct Particle
        instances.
        """
        instance = self.instances[test_name]
        base_particle = instance.base_particle
        nstates = instance.atomic_number + 1
        expected_particles = [Particle(base_particle, Z=Z) for Z in range(nstates)]
        assert expected_particles == instance.particles, (
            f"The expected Particle instances of {expected_particles} "
            f"are not all equal to the IonizationState particles of "
            f"{instance.particles} for test {test_name}."
        )

    def test_electron_density_from_n_elem_ionic_fractions(self):
        test_name = 'He'
        instance = self.instances[test_name]
        n_elem = test_cases[test_name]['n_elem']
        ionic_fractions = test_cases[test_name]['ionic_fractions']
        assert instance.n_elem == n_elem, \
            f"n_elem is not being stored correctly for test {test_name}"
        assert np.isclose(
            instance.n_e,
            np.sum(n_elem * ionic_fractions * np.array([0, 1, 2])),
            rtol=1e-12, atol=0 * u.m ** -3), \
            "n_e is not the expected value."


class Test_IonizationStates:

    @classmethod
    def setup_class(cls):

        cls.tests = {

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

        }

        cls.instances = {}

    tests = ['basic', 'quantities', 'just H', 'H acceptable error']

    @pytest.mark.parametrize('test', tests)
    def test_instantiation(self, test):
        try:
            self.instances[test] = IonizationStates(**self.tests[test])
        except Exception as exc:
            raise AtomicError(
                f"Cannot create IonizationStates instance for "
                f"test='{test}'") from exc

    @pytest.mark.parametrize('test', tests)
    def test_keys(self, test):
        input_keys = self.tests[test]['inputs'].keys()
        particles = [Particle(input_key) for input_key in input_keys]
        expected_keys = [p.particle for p in particles]
        actual_keys = [key for key in self.instances[test].ionic_fractions.keys()]
        elements = self.instances[test].elements

        assert actual_keys == expected_keys == elements, (
            f"For test='{test}', the following should be equal:\n"
            f"  actual_keys = {actual_keys}\n"
            f"expected_keys = {expected_keys}\n"
            f"     elements = {elements}"
        )

    @pytest.mark.parametrize('test', tests)
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

    @pytest.mark.parametrize('test', tests)
    def test_ionic_fractions(self, test):

        errmsg = ""

        elements = self.instances[test].elements
        input_keys = self.tests[test]["inputs"].keys()
        for element, input_key in zip(elements, input_keys):

            expected = np.array(self.tests[test]["inputs"][input_key])

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

    @pytest.mark.parametrize('test', tests)
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

    @pytest.mark.parametrize('test', tests)
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

    @pytest.mark.parametrize('test', tests)
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
