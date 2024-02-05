import collections
import itertools

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.particles import (
    atomic_number,
    atomic_symbol,
    charge_number,
    ionic_levels,
    isotope_symbol,
    particle_symbol,
)
from plasmapy.particles.exceptions import InvalidIsotopeError, ParticleError
from plasmapy.particles.ionization_state import IonicLevel, IonizationState
from plasmapy.particles.particle_class import Particle
from plasmapy.particles.particle_collections import ParticleList

ionic_fraction_table = [
    ("Fe 6+", 0.52, 5.2e-6 * u.m**-3),
    ("He 1+", None, None),
    ("H-2 0+", None, None),
]


@pytest.mark.parametrize(
    ("ion", "ionic_fraction", "number_density"), ionic_fraction_table
)
def test_ionic_level_attributes(ion, ionic_fraction, number_density) -> None:
    instance = IonicLevel(
        ion=ion, ionic_fraction=ionic_fraction, number_density=number_density
    )

    # Prepare to check for the default values when they are not set

    if ionic_fraction is None:
        ionic_fraction = np.nan
    if number_density is None:
        number_density = np.nan * u.m**-3

    assert Particle(ion) == Particle(instance.ionic_symbol)
    assert u.isclose(instance.ionic_fraction, ionic_fraction, equal_nan=True)
    assert u.isclose(instance.number_density, number_density, equal_nan=True)
    assert instance.charge_number == charge_number(ion)


@pytest.mark.parametrize(
    ("invalid_fraction", "expected_exception"),
    [(-1e-9, ParticleError), (1.00000000001, ParticleError), ("...", ParticleError)],
)
def test_ionic_level_invalid_inputs(invalid_fraction, expected_exception) -> None:
    """
    Test that IonicLevel raises exceptions when the ionic fraction
    is out of the interval [0,1] or otherwise invalid.
    """
    with pytest.raises(expected_exception):
        IonicLevel(ion="Fe 6+", ionic_fraction=invalid_fraction)


@pytest.mark.parametrize("invalid_particle", ["H", "e-", "Fe-56"])
def test_ionic_level_invalid_particles(invalid_particle) -> None:
    """
    Test that `~plasmapy.particles.IonicLevel` raises the appropriate
    exception when passed a particle that isn't a neutral or ion.
    """
    with pytest.raises(ParticleError):
        IonicLevel(invalid_particle, ionic_fraction=0)


@pytest.mark.parametrize(("ion1", "ion2"), [("Fe-56 6+", "Fe-56 5+"), ("H 1+", "D 1+")])
def test_ionic_level_comparison_with_different_ions(ion1, ion2) -> None:
    """
    Test that a `TypeError` is raised when an `IonicLevel` object
    is compared to an `IonicLevel` object of a different ion.
    """
    fraction = 0.251

    ionic_fraction_1 = IonicLevel(ion=ion1, ionic_fraction=fraction)
    ionic_fraction_2 = IonicLevel(ion=ion2, ionic_fraction=fraction)

    assert ionic_fraction_1 != ionic_fraction_2


def test_ionic_level_inequality_with_different_type() -> None:
    instance = IonicLevel("H 1+", 0.3)
    assert instance != "different type"


def test_ionization_state_ion_input_error() -> None:
    """
    Test that `~plasmapy.particles.IonizationState` raises the appropriate
    exception when an ion is the base particle and ionic fractions are
    specified
    """

    ion = "He 1+"
    unnecessary_ionic_fractions = [0.0, 0.0, 1.0]

    with pytest.raises(ParticleError):
        IonizationState(ion, ionic_fractions=unnecessary_ionic_fractions)


test_cases = {
    "Li": {
        "particle": "Li",
        "ionic_fractions": np.array([0.4, 0.3, 0.2, 0.1]),
        "tol": 1e-15,
    },
    "Li ground state": {
        "particle": "Li",
        "ionic_fractions": np.array([1, 0, 0, 0], dtype=int),
        "tol": 1e-15,
    },
    "H": {"particle": "H", "ionic_fractions": [0.6, 0.4], "tol": 1e-8},
    "H acceptable error": {
        "particle": "H",
        "ionic_fractions": [0.6, 0.400_000_001],
        "tol": 1e-8,
    },
    "D": {
        "particle": "deuterium",
        "ionic_fractions": [0.7, 0.3],
        "tol": 1e-15,
        "n_elem": 3e14 * u.m**-3,
    },
    "He": {
        "particle": "He",
        "ionic_fractions": [0.5, 0.3, 0.2],
        "n_elem": 1e20 * u.m**-3,
    },
    "number densities": {
        "particle": "T",
        "ionic_fractions": np.array([1e4, 1]) * u.cm**-3,
    },
}


@pytest.fixture(params=list(test_cases.values()), ids=list(test_cases.keys()))
def test_ionization_state(request):
    return IonizationState(**request.param)


def test_charge_numbers(test_ionization_state) -> None:
    """
    Test that an `IonizationState` instance has the correct charge
    numbers.
    """
    expected_charge_numbers = np.arange(test_ionization_state.atomic_number + 1)
    assert np.allclose(test_ionization_state.charge_numbers, expected_charge_numbers)


def test_equal_to_itself(He_ionization_state) -> None:
    """
    Test that `IonizationState.__eq__` returns `True for two identical
    `IonizationState` instances.
    """
    assert He_ionization_state == He_ionization_state  # noqa: PLR0124


@pytest.mark.parametrize(("tolerance", "output"), [(1e-8, True), (1e-9001, False)])
def test_equal_to_within_tolerance(tolerance: float, output) -> None:
    """
    Test that `IonizationState.__eq__` returns `True` for two
    `IonizationState` instances that differ within the inputted
    tolerance.
    """
    H = IonizationState(particle="H", ionic_fractions=[0.6, 0.4], tol=tolerance)
    H_acceptable_error = IonizationState(
        particle="H", ionic_fractions=[0.6, 0.400_000_001], tol=1e-8
    )
    assert (H == H_acceptable_error) == output


def test_equality_no_more_exception(test_ionization_state, He_ionization_state) -> None:
    """
    Test that comparisons of `IonizationState` instances for
    different elements does not fail.
    """
    if test_ionization_state == He_ionization_state:
        return
    assert test_ionization_state != He_ionization_state


def test_iteration(test_ionization_state) -> None:
    """Test that `IonizationState` instances iterate impeccably."""
    states = list(test_ionization_state)

    charge_numbers = [state.charge_number for state in states]
    ionic_symbols = [state.ionic_symbol for state in states]

    try:
        base_symbol = isotope_symbol(ionic_symbols[0])
    except InvalidIsotopeError:
        base_symbol = atomic_symbol(ionic_symbols[0])
    finally:
        atomic_numb = atomic_number(ionic_symbols[1])

    errors = []

    expected_charges = np.arange(atomic_numb + 1)
    if not np.all(charge_numbers == expected_charges):
        errors.append(
            f"The resulting charge numbers are {charge_numbers}, "
            f"which are not equal to the expected charge numbers, "
            f"which are {expected_charges}."
        )

    np.testing.assert_allclose(sum(test_ionization_state.ionic_fractions), 1)

    expected_particles = [Particle(base_symbol, Z=charge) for charge in charge_numbers]
    expected_symbols = [particle.ionic_symbol for particle in expected_particles]
    if ionic_symbols != expected_symbols:
        errors.append(
            f"The resulting ionic symbols are {ionic_symbols}, "
            f"which are not equal to the expected ionic symbols of "
            f"{expected_symbols}."
        )

    if errors:
        errors.insert(
            0,
            (
                f"The test of {test_ionization_state} has "
                f"resulted in the following errors when attempting to "
                f"iterate."
            ),
        )
        errmsg = " ".join(errors)
        pytest.fail(errmsg)


def test_normalization() -> None:
    """
    Test that `_is_normalized` returns `False` when there is an
    error greater than the tolerance, and `True` after normalizing.
    """
    H = IonizationState(particle="H", ionic_fractions=[0.6, 0.400_000_001], tol=1e-8)
    assert not H._is_normalized(tol=1e-15)
    H.normalize()
    assert H._is_normalized(tol=1e-15)


def test_identifications(test_ionization_state) -> None:
    """
    Test that the identification attributes for test
    `IonizationState` instances match the expected values from the
    `Particle` instance.
    """

    Identifications = collections.namedtuple(
        "Identifications", ["element", "isotope", "atomic_number"]
    )

    expected_identifications = Identifications(
        test_ionization_state.element,
        test_ionization_state.isotope,
        test_ionization_state.atomic_number,
    )

    expected_element = test_ionization_state._particle.element
    expected_isotope = test_ionization_state._particle.isotope
    expected_atomic_number = test_ionization_state._particle.atomic_number

    resulting_identifications = Identifications(
        expected_element, expected_isotope, expected_atomic_number
    )

    assert resulting_identifications == expected_identifications, (
        f"For {test_ionization_state}, the resulting "
        f"identifications of {resulting_identifications} differ "
        f"from the expected identifications of "
        f"{expected_identifications}."
    )


def test_as_particle_list(test_ionization_state) -> None:
    """
    Test that `IonizationState` returns the correct `Particle`
    instances.
    """
    atom = test_ionization_state.base_particle
    nstates = test_ionization_state.atomic_number + 1
    expected_particles = [Particle(atom, Z=Z) for Z in range(nstates)]
    actual_particles = test_ionization_state.to_list()
    assert expected_particles == actual_particles


def test_getitem(test_ionization_state) -> None:
    """
    Test that `IonizationState.__getitem__` returns the same value
    when using equivalent keys (charge number, particle symbol, and
    `Particle` instance).

    For example, if we create

    >>> He_states = IonizationState("He", [0.2, 0.3, 0.5])

    then this checks to make sure that `He_states[2]`,
    `He_states['He 2+']`, and `He_states[Particle('He 2+')]` all
    return the same result.

    """
    particle_name = test_ionization_state.base_particle

    charge_numbers = np.arange(test_ionization_state.atomic_number + 1)
    symbols = [particle_symbol(particle_name, Z=Z) for Z in charge_numbers]
    particles = test_ionization_state.to_list()

    errors = []

    # In the following loop, instance[key] will return a namedtuple
    # or class which may contain Quantity objects with values of
    # numpy.nan.  Because of the difficulty of comparing nans in
    # these objects, we compare the string representations instead
    # (see Astropy issue #7901 on GitHub).

    for keys in zip(charge_numbers, symbols, particles):
        set_of_str_values = {str(test_ionization_state[key]) for key in keys}
        if len(set_of_str_values) != 1:
            errors.append(
                f"\n\n"
                f"The following keys for {test_ionization_state} did not "
                f"produce identical outputs as required: {keys}. "
                f"The set containing string representations of"
                f"the values is:\n\n{set_of_str_values}"
            )

    if errors:
        pytest.fail(str.join("", errors))


@pytest.fixture()
def He_ionization_state():
    return IonizationState(
        particle="He",
        ionic_fractions=[0.5, 0.3, 0.2],
        n_elem=1e20 * u.m**-3,
    )


@pytest.mark.parametrize("index", [-1, 4, "Li"])
def test_indexing_error(He_ionization_state, index) -> None:
    """
    Test that an `IonizationState` instance cannot be indexed
    outside of the bounds of allowed charge numbers.
    """
    with pytest.raises(ParticleError):
        He_ionization_state[index]


def test_inequality_with_different_type(He_ionization_state) -> None:
    assert He_ionization_state != "different type"


@pytest.mark.parametrize("tol", [-1e-16, 1.0000001])
def test_invalid_tolerances(He_ionization_state, tol: float) -> None:
    """Test that invalid tolerances raise appropriate errors."""
    with pytest.raises(ValueError):
        He_ionization_state.tol = tol


def test_State_equality_and_getitem(He_ionization_state) -> None:
    charge = 2
    symbol = "He 2+"
    result_from_charge = He_ionization_state[charge]
    result_from_symbol = He_ionization_state[symbol]
    assert result_from_charge == result_from_symbol


ions = ["Fe 6+", "p", "He-4 0+", "triton", "alpha", "Ne +0"]


@pytest.mark.parametrize("ion", ions)
def test_IonizationState_ionfracs_from_ion_input(ion) -> None:
    ionization_state = IonizationState(ion)
    ion_particle = Particle(ion)
    actual_ionic_fractions = ionization_state.ionic_fractions

    expected_ionic_fractions = np.zeros(ion_particle.atomic_number + 1)
    expected_ionic_fractions[ion_particle.charge_number] = 1.0

    np.testing.assert_allclose(
        expected_ionic_fractions,
        actual_ionic_fractions,
        atol=1e-16,
        err_msg=f"The returned ionic fraction for IonizationState({ion!r}) "
        f"should have entirely been in the Z = {ion_particle.charge_number} "
        f"level.",
    )


@pytest.mark.parametrize("ion", ions)
def test_IonizationState_base_particles_from_ion_input(ion) -> None:
    """
    Test that supplying an ion to IonizationState will result in the
    base particle being the corresponding isotope or ion and that the
    ionic fraction of the corresponding charge level is 100%.
    """

    ionization_state = IonizationState(ion)
    ion_particle = Particle(ion)

    expected_base_particle = ion_particle.isotope or ion_particle.element
    if expected_base_particle != ionization_state.base_particle:
        pytest.fail(
            f"The expected base particle was {expected_base_particle}, "
            f"but the returned base particle was {ionization_state.base_particle}. "
        )


expected_properties = {
    "T_e": 5000.0 * u.K,
    "isotope": "He-4",
    "element": "He",
    "atomic_number": 2,
    "Z_mean": 1.3,
    "Z_rms": 1.51657508881031,
    "n_e": 1.3e19 * u.m**-3,
    "n_elem": 1e19 * u.m**-3,
    "charge_numbers": [0, 1, 2],
    "ionic_fractions": np.array([0.2, 0.3, 0.5]),
    "ionic_symbols": ["He-4 0+", "He-4 1+", "He-4 2+"],
    "number_densities": np.array([2e18, 3e18, 5e18]) * u.m**-3,
    "tol": 2e-14,
}


@pytest.fixture()
def instance():
    kwargs = {
        "particle": "He-4",
        "ionic_fractions": [0.2, 0.3, 0.5],
        "T_e": 5.0 * u.kK,
        "tol": 2e-14,
        "n_elem": 1e13 * u.cm**-3,
    }

    return IonizationState(**kwargs)


@pytest.mark.parametrize("key", expected_properties)
def test_IonizationState_attributes(instance, key) -> None:
    """
    Test a specific case that the `IonizationState` attributes are
    working as expected.
    """
    expected = expected_properties[key]
    actual = getattr(instance, key)

    if isinstance(expected, u.Quantity):
        assert expected.unit == actual.unit, f"Unit mismatch for IonizationState.{key}"
        assert np.allclose(
            expected, actual, atol=1e-15 * expected.unit
        ), f"Quantity.value mismatch for IonizationState.{key}"
    else:
        try:
            assert expected == actual
        except ValueError:
            assert np.allclose(expected, actual)


def test_IonizationState_methods(instance) -> None:
    assert instance._is_normalized()
    assert str(instance) == "<IonizationState instance for He-4>"


def test_IonizationState_ion_temperatures(instance) -> None:
    for ionic_level in instance:
        assert instance.T_e == ionic_level.T_i


@pytest.mark.xfail(
    reason="IonizationState currently does not store IonicLevels, but generates them on the fly!"
)
def test_IonizationState_ion_temperature_persistence(instance) -> None:
    instance[0].T_i += 1 * u.K
    assert instance[0].T_i - instance.T_e == (1 * u.K)


@pytest.mark.parametrize(
    "T_i",
    [
        10 * u.eV,
        1000 * u.K,
        None,
        u.Quantity([1, 1, 10], u.eV),
        u.Quantity([1000, 1000, 10000], u.K),
    ],
)
def test_set_T_i(instance, T_i) -> None:
    instance.T_i = T_i


def test_default_T_i_is_T_e(instance) -> None:
    T_i = instance.T_i
    assert_quantity_allclose(T_i, instance.T_e)
    assert len(T_i) == 3


@pytest.mark.parametrize(
    ("T_i", "expectation"),
    [
        (10 * u.m, pytest.raises(u.UnitTypeError)),
        (
            u.Quantity([1, 1], u.eV),
            pytest.raises(
                ParticleError,
                match="common temperature for all ions, or a set of 3 of them",
            ),
        ),
        (
            u.Quantity([1] * 5, u.eV),
            pytest.raises(ParticleError, match="five is right out"),
        ),
    ],
)
def test_set_T_i_with_errors(instance, T_i, expectation) -> None:
    with expectation:
        instance.T_i = T_i


def test_slicing(instance) -> None:
    instance[1:]


def test_len(instance) -> None:
    assert len(instance) == 3


def test_nans() -> None:
    """
    Test that when no ionic fractions or temperature are inputted,
    the result is an array full of `~numpy.nan` of the right size.
    """
    element = "He"
    nstates = atomic_number(element) + 1
    instance = IonizationState(element)
    assert (
        len(instance.ionic_fractions) == nstates
    ), f"Incorrect number of ionization states for {element}"
    assert np.all([np.isnan(instance.ionic_fractions)]), (
        "The ionic fractions for IonizationState are not defaulting "
        "to numpy.nan when not set by user."
    )


def test_setting_ionic_fractions() -> None:
    """
    Test the setter for the ``ionic_fractions`` attribute on
    `~plasmapy.particles.IonizationState`.
    """
    instance = IonizationState("He")
    new_ionic_fractions = [0.2, 0.5, 0.3]
    instance.ionic_fractions = new_ionic_fractions
    assert np.allclose(instance.ionic_fractions, new_ionic_fractions)


class Test_IonizationStateNumberDensitiesSetter:
    """Test that setting IonizationState.number_densities works correctly."""

    def setup_class(self) -> None:
        self.element = "H"
        self.valid_number_densities = u.Quantity([0.1, 0.2], unit=u.m**-3)
        self.expected_n_elem = np.sum(self.valid_number_densities)
        self.expected_ionic_fractions = (
            self.valid_number_densities / self.expected_n_elem
        )
        try:
            self.instance = IonizationState(self.element)
        except Exception:  # noqa: BLE001
            pytest.fail(
                "Unable to instantiate IonizationState with no ionic fractions."
            )

    def test_setting_number_densities(self) -> None:
        try:
            self.instance.number_densities = self.valid_number_densities
        except Exception:  # noqa: BLE001
            pytest.fail(
                f"Unable to set number densities of {self.element} to "
                f"{self.valid_number_densities}."
            )

        assert u.quantity.allclose(
            self.instance.number_densities, self.valid_number_densities
        ), (
            f"The number densities of {self.element} were set to "
            f"{self.instance.number_densities} instead of the expected "
            f"value of {self.valid_number_densities}."
        )

    def test_ionic_fractions(self) -> None:
        assert np.allclose(
            self.instance.ionic_fractions, self.expected_ionic_fractions
        ), (
            "The IonizationState.ionic_fractions attribute was not set "
            "correctly after the number densities were set."
        )

    def test_n_elem(self) -> None:
        assert u.quantity.allclose(self.instance.n_elem, self.expected_n_elem), (
            "IonizationState.n_elem not set correctly after "
            "number_densities was set."
        )

    def test_n_e(self) -> None:
        assert u.quantity.allclose(
            self.instance.n_e, self.valid_number_densities[1]
        ), "IonizationState.n_e not set correctly after number_densities was set."

    def test_that_negative_density_raises_error(self) -> None:
        with pytest.raises(ParticleError, match="cannot be negative"):
            self.instance.number_densities = u.Quantity([-0.1, 0.2], unit=u.m**-3)

    def test_incorrect_number_of_charge_states_error(self) -> None:
        with pytest.raises(ParticleError, match="Incorrect number of charge states"):
            self.instance.number_densities = u.Quantity([0.1, 0.2, 0.3], unit=u.m**-3)

    def test_incorrect_units_error(self) -> None:
        with pytest.raises(u.UnitsError):
            self.instance.number_densities = u.Quantity([0.1, 0.2], unit=u.kg)

    # The following two tests are not related to setting the
    # number_densities attribute, but are helpful to test anyway.

    def test_T_e_isnan_when_not_set(self) -> None:
        assert np.isnan(self.instance.T_e)

    def test_kappa_isinf_when_not_set(self) -> None:
        assert np.isinf(self.instance.kappa)


def test_iteration_with_nested_iterator() -> None:
    hydrogen = IonizationState("p+", n_elem=1e20 * u.m**-3, T_e=10 * u.eV)

    i = sum(1 for _, __ in itertools.product(hydrogen, hydrogen))
    assert i == 4


def test_ionization_state_inequality_and_identity() -> None:
    deuterium_states = IonizationState("D+", n_elem=1e20 * u.m**-3, T_e=10 * u.eV)
    tritium_states = IonizationState("T+", n_elem=1e20 * u.m**-3, T_e=10 * u.eV)
    assert deuterium_states != tritium_states


physical_properties = ["charge", "mass"]

particles_and_ionfracs = [
    ("H-1", np.array([1, 0])),
    ("H-1", np.array([0, 1])),
    ("H-1", np.array([0.3, 0.7])),
    ("He-4", np.array([0.2, 0.5, 0.3])),
    ("Li-7", np.array([0.21, 0.01, 0.28, 0.5])),
]


@pytest.mark.parametrize("physical_property", physical_properties)
@pytest.mark.parametrize(("base_particle", "ionic_fractions"), particles_and_ionfracs)
def test_weighted_mean_ion(base_particle, ionic_fractions, physical_property) -> None:
    """
    Test that `IonizationState.average_ion` gives a |CustomParticle|
    instance with the expected mass or charge when calculating the
    weighted mean.
    """
    ionization_state = IonizationState(base_particle, ionic_fractions)
    ions = ionic_levels(base_particle)
    physical_quantity = getattr(ions, physical_property)
    expected_mean_quantity = np.average(physical_quantity, weights=ionic_fractions)
    mean_ion = ionization_state.average_ion()
    actual_mean_quantity = getattr(mean_ion, physical_property)
    assert_quantity_allclose(actual_mean_quantity, expected_mean_quantity)


@pytest.mark.parametrize("physical_property", physical_properties)
@pytest.mark.parametrize(("base_particle", "ionic_fractions"), particles_and_ionfracs)
def test_weighted_rms_ion(base_particle, ionic_fractions, physical_property) -> None:
    """
    Test that `IonizationState.average_ion` gives a |CustomParticle|
    instances with the expected mass or charge when calculating the
    weighted root mean square.
    """
    ionization_state = IonizationState(base_particle, ionic_fractions)
    ions = ionic_levels(base_particle)
    physical_quantity = getattr(ions, physical_property)
    expected_rms_quantity = np.sqrt(
        np.average(physical_quantity**2, weights=ionic_fractions)
    )
    kwargs = {f"use_rms_{physical_property}": True}
    rms_ion = ionization_state.average_ion(**kwargs)
    actual_rms_quantity = getattr(rms_ion, physical_property)
    assert_quantity_allclose(actual_rms_quantity, expected_rms_quantity)


def test_exclude_neutrals_from_average_ion() -> None:
    """
    Test that the `IonizationState.average_ion` method returns a
    |CustomParticle| that does not include neutrals in the averaging
    when the ``include_neutrals`` keyword is `False`.
    """
    base_particle = Particle("He-4")
    ionization_state_without_neutrals = IonizationState(base_particle, [0, 0.2, 0.8])
    expected_average_ion = ionization_state_without_neutrals.average_ion()
    ionization_state_with_neutrals = IonizationState(base_particle, [0.50, 0.1, 0.4])
    actual_average_ion = ionization_state_with_neutrals.average_ion(
        include_neutrals=False
    )
    assert actual_average_ion == expected_average_ion


@pytest.mark.parametrize("physical_property", physical_properties)
def test_comparison_to_equivalent_particle_list(physical_property) -> None:
    """
    Test that `IonizationState.average_ion` gives consistent results with
    `ParticleList.average_particle` when the ratios of different particles
    is the same between the `IonizationState` and the `ParticleList`.
    """
    particles = ParticleList(2 * ["He-4 0+"] + 3 * ["He-4 1+"] + 5 * ["He-4 2+"])
    ionization_state = IonizationState("He-4", [0.2, 0.3, 0.5])
    kwargs = {f"use_rms_{physical_property}": True}
    expected_average_particle = particles.average_particle(**kwargs)
    expected_average_quantity = getattr(expected_average_particle, physical_property)
    actual_average_particle = ionization_state.average_ion(**kwargs)
    actual_average_quantity = getattr(actual_average_particle, physical_property)
    assert_quantity_allclose(actual_average_quantity, expected_average_quantity)
