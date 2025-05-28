import collections
import inspect
import io
import json

import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest
from astropy.constants import c, e, m_e, m_n, m_p

from plasmapy.particles import json_load_particle, json_loads_particle, molecule
from plasmapy.particles._isotopes import data_about_isotopes
from plasmapy.particles.atomic import known_isotopes
from plasmapy.particles.exceptions import (
    ChargeError,
    InvalidElementError,
    InvalidIonError,
    InvalidIsotopeError,
    InvalidParticleError,
    MissingParticleDataError,
    MissingParticleDataWarning,
    ParticleError,
    ParticleWarning,
)
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
    ParticleLike,
    valid_categories,
)
from plasmapy.utils import roman
from plasmapy.utils._pytest_helpers import run_test_equivalent_calls
from plasmapy.utils.code_repr import call_string

# (arg, kwargs, results_dict)
test_Particle_table = [
    (
        "neutron",
        {},
        {
            "symbol": "n",
            "element": None,
            "isotope": None,
            "isotope_name": InvalidElementError,
            "ionic_symbol": None,
            "roman_symbol": None,
            "is_ion": False,
            "is_electron": False,
            "charge_number": 0,
            "atomic_number": InvalidElementError,
            "mass_number": InvalidIsotopeError,
            "baryon_number": 1,
            "lepton_number": 0,
            "mass": m_n,
            "nuclide_mass": m_n,
            "nuclear_binding_energy": 0 * u.J,
            "periodic_table.group": InvalidElementError,
        },
    ),
    (
        "p+",
        {},
        {
            "symbol": "p+",
            "element": "H",
            "element_name": "hydrogen",
            "isotope": "H-1",
            "isotope_name": "hydrogen-1",
            "ionic_symbol": "p+",
            "roman_symbol": "H-1 II",
            "is_ion": True,
            "mass": m_p,
            "nuclide_mass": m_p,
            "nucleus": Particle("p+"),
            "charge_number": 1,
            "charge.value": e.si.value,
            "spin": 1 / 2,
            "half_life": np.inf * u.s,
            "atomic_number": 1,
            "mass_number": 1,
            "lepton_number": 0,
            "baryon_number": 1,
            "__str__()": "p+",
            "__repr__()": 'Particle("p+")',
            'is_category("fermion")': True,
            'is_category(["fermion"])': True,
            'is_category({"fermion"})': True,
            'is_category(any_of=("boson", "fermion"))': True,
            'is_category(require=["boson", "fermion"])': False,
            'is_category(("element", "isotope", "ion"))': True,
            'is_category("charged")': True,
            "periodic_table.group": 1,
            "periodic_table.block": "s",
            "periodic_table.period": 1,
            "periodic_table.category": "nonmetal",
            "nuclear_binding_energy": 0 * u.J,
            "recombine()": "H-1 0+",
        },
    ),
    (
        "H",  # proton
        {"Z": 1, "mass_numb": 1},
        {
            "symbol": "p+",
            "element": "H",
            "element_name": "hydrogen",
            "isotope": "H-1",
            "isotope_name": "hydrogen-1",
            "ionic_symbol": "p+",
            "roman_symbol": "H-1 II",
            "is_ion": True,
            "mass": m_p,
            "nuclide_mass": m_p,
            "nucleus": Particle("p+"),
            "charge_number": 1,
            "charge.value": e.si.value,
            "spin": 1 / 2,
            "half_life": np.inf * u.s,
            "atomic_number": 1,
            "mass_number": 1,
            "lepton_number": 0,
            "baryon_number": 1,
            "__str__()": "p+",
            "__repr__()": 'Particle("p+")',
            'is_category("fermion")': True,
            'is_category(["fermion"])': True,
            'is_category({"fermion"})': True,
            'is_category(any_of=("boson", "fermion"))': True,
            'is_category(require=["boson", "fermion"])': False,
            'is_category(("element", "isotope", "ion"))': True,
            'is_category("charged")': True,
            "periodic_table.group": 1,
            "periodic_table.block": "s",
            "periodic_table.period": 1,
            "periodic_table.category": "nonmetal",
            "nuclear_binding_energy": 0 * u.J,
            "recombine()": "H-1 0+",
        },
    ),
    (
        "H",  # without charge or mass number information
        {},
        {
            "symbol": "H",
            "element": "H",
            "isotope": None,
            "isotope_name": InvalidIsotopeError,
            "ionic_symbol": None,
            "roman_symbol": ChargeError,
            "is_ion": False,
            "charge": np.nan * u.C,
            "charge_number": ChargeError,
            "mass_number": InvalidIsotopeError,
            "baryon_number": ParticleError,
            "lepton_number": 0,
            "half_life": InvalidIsotopeError,
            "standard_atomic_weight": (1.008 * u.u).to(u.kg),
            "mass": (1.008 * u.u).to(u.kg),
            "nuclide_mass": InvalidIsotopeError,
            'is_category("charged")': False,
            'is_category("nonmetal")': True,
            'is_category("proton")': False,
        },
    ),
    (
        "p-",
        {},
        {
            "symbol": "p-",
            "element": None,
            "element_name": InvalidElementError,
            "isotope": None,
            "isotope_name": InvalidElementError,
            "ionic_symbol": None,
            "roman_symbol": None,
            "is_ion": False,
            "mass": m_p,
            "charge_number": -1,
            "spin": 1 / 2,
            "half_life": np.inf * u.s,
            "atomic_number": InvalidElementError,
            "mass_number": InvalidIsotopeError,
            "lepton_number": 0,
            "baryon_number": -1,
            "__str__()": "p-",
            "__repr__()": 'Particle("p-")',
            "periodic_table.group": InvalidElementError,
        },
    ),
    (
        "e-",
        {},
        {
            "symbol": "e-",
            "element": None,
            "element_name": InvalidElementError,
            "isotope": None,
            "isotope_name": InvalidElementError,
            "ionic_symbol": None,
            "roman_symbol": None,
            "is_ion": False,
            "mass": m_e,
            "charge_number": -1,
            "spin": 1 / 2,
            "half_life": np.inf * u.s,
            "atomic_number": InvalidElementError,
            "lepton_number": 1,
            "baryon_number": 0,
            "__str__()": "e-",
            "__repr__()": 'Particle("e-")',
            "nuclear_binding_energy": InvalidIsotopeError,
            "periodic_table.group": InvalidElementError,
            "periodic_table.block": InvalidElementError,
            "periodic_table.period": InvalidElementError,
            "periodic_table.category": InvalidElementError,
        },
    ),
    (
        "e+",
        {},
        {
            "symbol": "e+",
            "element": None,
            "isotope": None,
            "isotope_name": InvalidElementError,
            "ionic_symbol": None,
            "roman_symbol": None,
            "is_ion": False,
            "mass": m_e,
            "mass_energy": (m_e * c**2).to("J"),
            "nuclide_mass": InvalidIsotopeError,
            "charge_number": 1,
            "spin": 1 / 2,
            "half_life": np.inf * u.s,
            "atomic_number": InvalidElementError,
            "lepton_number": -1,
            "baryon_number": 0,
            'is_category(require="positron")': True,
            'is_category(any_of={"positron"})': True,
            'is_category(exclude="positron")': False,
            'is_category("ion")': False,
            'is_category("invalid_category")': ParticleError,
            "__str__()": "e+",
            "__repr__()": 'Particle("e+")',
            "periodic_table.group": InvalidElementError,
            "periodic_table.block": InvalidElementError,
            "periodic_table.period": InvalidElementError,
            "periodic_table.category": InvalidElementError,
        },
    ),
    (
        "H 1-",
        {},
        {
            "symbol": "H 1-",
            "element": "H",
            "isotope": None,
            "isotope_name": InvalidIsotopeError,
            "ionic_symbol": "H 1-",
            "roman_symbol": roman.OutOfRangeError,
            "is_ion": True,
            "charge_number": -1,
            "mass_number": InvalidIsotopeError,
            "baryon_number": MissingParticleDataError,
            "lepton_number": 0,
            "half_life": InvalidIsotopeError,
            "standard_atomic_weight": InvalidElementError,
            "nuclide_mass": InvalidIsotopeError,
            'is_category("charged")': True,
            'is_category("nonmetal")': True,
            'is_category("proton")': False,
        },
    ),
    (
        "H-1 0+",
        {},
        {
            "symbol": "H-1 0+",
            "element": "H",
            "isotope": "H-1",
            "isotope_name": "hydrogen-1",
            "ionic_symbol": "H-1 0+",
            "roman_symbol": "H-1 I",
            "is_ion": False,
            "charge": 0 * u.C,
            "charge_number": 0,
            "mass_number": 1,
            "baryon_number": 1,
            "lepton_number": 0,
            "half_life": np.inf * u.s,
            "nuclide_mass": m_p,
            "nucleus": Particle("p+"),
            'is_category("charged")': False,
            'is_category("uncharged")': True,
            'is_category("ion")': False,
            'is_category("nonmetal")': True,
            'is_category("proton")': False,
        },
    ),
    (
        "D+",
        {},
        {
            "symbol": "D 1+",
            "element": "H",
            "element_name": "hydrogen",
            "isotope": "D",
            "isotope_name": "deuterium",
            "ionic_symbol": "D 1+",
            "roman_symbol": "D II",
            "is_ion": True,
            "charge_number": 1,
            "atomic_number": 1,
            "mass_number": 2,
            "baryon_number": 2,
            "lepton_number": 0,
            "neutron_number": 1,
            "nucleus": Particle("D 1+"),
            'is_category(require=("ion", "isotope"))': True,
            "periodic_table.group": 1,
            "periodic_table.block": "s",
            "periodic_table.period": 1,
            "periodic_table.category": "nonmetal",
        },
    ),
    (
        "tritium",
        {"Z": 1},
        {
            "symbol": "T 1+",
            "element": "H",
            "isotope": "T",
            "isotope_name": "tritium",
            "ionic_symbol": "T 1+",
            "roman_symbol": "T II",
            "is_ion": True,
            "charge_number": 1,
            "atomic_number": 1,
            "mass_number": 3,
            "baryon_number": 3,
            "lepton_number": 0,
            "neutron_number": 2,
            "nucleus": Particle("T 1+"),
            'is_category("ion", "isotope")': True,
            'is_category(require="uncharged")': False,
            "periodic_table.group": 1,
        },
    ),
    (
        "Fe",
        {"Z": 17, "mass_numb": 56},
        {
            "symbol": "Fe-56 17+",
            "element": "Fe",
            "element_name": "iron",
            "isotope": "Fe-56",
            "isotope_name": "iron-56",
            "ionic_symbol": "Fe-56 17+",
            "roman_symbol": "Fe-56 XVIII",
            "is_electron": False,
            "is_ion": True,
            "charge_number": 17,
            "atomic_number": 26,
            "mass_number": 56,
            "baryon_number": 56,
            "__str__()": "Fe-56 17+",
            "__repr__()": 'Particle("Fe-56 17+")',
            'is_category("element")': True,
            'is_category("ion")': True,
            'is_category("isotope")': True,
        },
    ),
    (
        "alpha",
        {},
        {
            "symbol": "He-4 2+",
            "element": "He",
            "element_name": "helium",
            "isotope": "He-4",
            "isotope_name": "helium-4",
            "ionic_symbol": "He-4 2+",
            "roman_symbol": "He-4 III",
            "mass_energy": 5.971919969131517e-10 * u.J,
            "is_ion": True,
            "charge_number": 2,
            "atomic_number": 2,
            "mass_number": 4,
            "baryon_number": 4,
            "lepton_number": 0,
            "nucleus": Particle("He-4 2+"),
            "half_life": np.inf * u.s,
            "recombine()": Particle("He-4 1+"),
        },
    ),
    (
        "He-4 0+",
        {},
        {
            "symbol": "He-4 0+",
            "element": "He",
            "isotope": "He-4",
            "nucleus": Particle("He-4 2+"),
            "mass_energy": 5.971919969131517e-10 * u.J,
        },
    ),
    (
        "Li",
        {"mass_numb": 7},
        {
            "symbol": "Li-7",
            "element": "Li",
            "element_name": "lithium",
            "isotope": "Li-7",
            "isotope_name": "lithium-7",
            "ionic_symbol": None,
            "roman_symbol": ChargeError,
            "is_ion": False,
            "charge_number": ChargeError,
            "atomic_number": 3,
            "mass_number": 7,
            "neutron_number": 4,
            "baryon_number": 7,
            "half_life": np.inf * u.s,
            "nuclide_mass": 1.1647614796180463e-26 * u.kg,
            "nucleus": Particle("Li-7 3+"),
        },
    ),
    (
        "Cn-276",
        {"Z": 22},
        {
            "symbol": "Cn-276 22+",
            "element": "Cn",
            "isotope": "Cn-276",
            "isotope_name": "copernicium-276",
            "ionic_symbol": "Cn-276 22+",
            "roman_symbol": "Cn-276 XXIII",
            "is_ion": True,
            "element_name": "copernicium",
            "charge_number": 22,
            "atomic_number": 112,
            "mass_number": 276,
            "neutron_number": 164,
            "baryon_number": 276,
            "lepton_number": 0,
        },
    ),
    (
        "muon",
        {},
        {
            "symbol": "mu-",
            "element": None,
            "isotope": None,
            "isotope_name": InvalidElementError,
            "ionic_symbol": None,
            "roman_symbol": None,
            "is_ion": False,
            "charge_number": -1,
            "atomic_number": InvalidElementError,
            "mass_number": InvalidIsotopeError,
            "baryon_number": 0,
            "lepton_number": 1,
        },
    ),
    (
        "nu_tau",
        {},
        {
            "symbol": "nu_tau",
            "element": None,
            "isotope": None,
            "isotope_name": InvalidElementError,
            "mass": np.nan * u.kg,
            "mass_energy": np.nan * u.J,
            "charge_number": 0,
            "mass_number": InvalidIsotopeError,
            "element_name": InvalidElementError,
            "baryon_number": 0,
            "lepton_number": 1,
            "half_life": np.inf * u.s,
            "is_electron": False,
            "is_ion": False,
            'is_category("fermion")': True,
            'is_category("neutrino")': True,
            'is_category("boson")': False,
            'is_category("matter", exclude={"antimatter"})': True,
            'is_category("matter", exclude=["antimatter"])': True,
            'is_category("matter", exclude="antimatter")': True,
            'is_category(any_of={"matter", "boson"})': True,
            'is_category(any_of=["antimatter", "boson", "charged"])': False,
            'is_category(["fermion", "lepton"], exclude="matter")': False,
            'is_category("lepton", "invalid")': ParticleError,
            'is_category(["boson"], exclude=["lepton", "invalid"])': ParticleError,
            'is_category("boson", exclude="boson")': ParticleError,
            'is_category(any_of="boson", exclude="boson")': ParticleError,
        },
    ),
    (Particle("C"), {}, {"symbol": "C", "atomic_number": 6, "element": "C"}),
    (
        Particle("C"),
        {"Z": 3, "mass_numb": 14},
        {
            "symbol": "C-14 3+",
            "element": "C",
            "isotope": "C-14",
            "ionic_symbol": "C-14 3+",
        },
    ),
]


@pytest.mark.parametrize(("arg", "kwargs", "expected_dict"), test_Particle_table)
def test_Particle_class(arg, kwargs, expected_dict):
    """
    Test that `~plasmapy.particles.Particle` objects for different
    subatomic particles, elements, isotopes, and ions return the
    expected properties.  Provide a detailed error message that lists
    all of the inconsistencies with the expected results.
    """

    call = call_string(Particle, arg, kwargs)
    errmsg = ""

    try:
        particle = Particle(arg, **kwargs)  # noqa: F841
    except Exception as exc:
        raise ParticleError(f"Problem creating {call}") from exc

    for key in expected_dict:
        expected = expected_dict[key]

        if inspect.isclass(expected) and issubclass(expected, Exception):
            # Exceptions are expected to be raised when accessing certain
            # attributes for some particles.  For example, accessing a
            # neutrino's mass should raise a MissingParticleDataError since
            # only upper limits of neutrino masses are presently available.
            # If expected_dict[key] is an exception, then check to make
            # sure that this exception is raised.

            try:
                with pytest.raises(expected):
                    exec(f"particle.{key}")  # noqa: S102
            except pytest.fail.Exception:
                errmsg += f"\n{call}[{key}] does not raise {expected}."
            except Exception:  # noqa: BLE001
                errmsg += (
                    f"\n{call}[{key}] does not raise {expected} but "
                    f"raises a different exception."
                )

        else:
            try:
                result = eval(f"particle.{key}")  # noqa: S307
                assert result == expected or u.isclose(result, expected, equal_nan=True)
            except AssertionError:
                errmsg += (
                    f"\n{call}.{key} returns {result} instead "
                    f"of the expected value of {expected}."
                )
            except Exception:  # noqa: BLE001
                errmsg += f"\n{call}.{key} raises an unexpected exception."

    if errmsg:
        raise Exception(f"Problems with {call}:{errmsg}")  # noqa: TRY002


equivalent_particles_table: list[list[ParticleLike]] = [
    ["H", "hydrogen", "hYdRoGeN"],
    ["p+", "proton", "H-1+", "H-1 1+", "H-1 +1"],
    ["D", "H-2", "Hydrogen-2", "deuterium"],
    ["T", "H-3", "Hydrogen-3", "tritium"],
    ["alpha", "He-4++", "He-4 2+", "He-4 +2", "He-4 III"],
    ["e-", "electron", "e"],
    ["e+", "positron"],
    ["p-", "antiproton"],
    ["n", "n-1", "neutron", "NEUTRON"],
    ["muon", "mu-", "muon-"],
    ["tau", "tau-"],
    [Particle("Fe 5+"), Particle("Fe 4+").ionize()],
    [Particle("He-4 0+"), Particle("alpha").recombine(2)],
]


@pytest.mark.parametrize("equivalent_particles", equivalent_particles_table)
def test_Particle_equivalent_cases(equivalent_particles) -> None:
    """Test that all instances of a list of particles are equivalent."""
    run_test_equivalent_calls(Particle, *equivalent_particles)


# args, kwargs, attribute, exception
test_Particle_error_table = [
    (["a"], {}, "", InvalidParticleError),
    (["d+"], {"mass_numb": 9}, "", InvalidParticleError),
    (["H"], {"mass_numb": 99}, "", InvalidParticleError),
    (["Au-818"], {}, "", InvalidParticleError),
    (["Au-12"], {}, "", InvalidParticleError),
    (["Au"], {"mass_numb": 13}, "", InvalidParticleError),
    (["Au"], {"mass_numb": 921}, "", InvalidParticleError),
    (["e-"], {"Z": -1}, "", InvalidParticleError),
    (["e-"], {}, ".atomic_number", InvalidElementError),
    (["alpha"], {}, ".standard_atomic_weight", InvalidElementError),
    (["Fe-56"], {}, ".standard_atomic_weight", InvalidElementError),
    (["e-"], {}, ".standard_atomic_weight", InvalidElementError),
    (["tau-"], {}, ".element_name", InvalidElementError),
    (["tau+"], {}, ".atomic_number", InvalidElementError),
    (["neutron"], {}, ".atomic_number", InvalidElementError),
    (["H"], {"Z": 0}, ".mass_number", InvalidIsotopeError),
    (["neutron"], {}, ".mass_number", InvalidIsotopeError),
    (["He"], {"mass_numb": 4}, ".charge_number", ChargeError),
    (["Fe"], {}, ".spin", MissingParticleDataError),
    ([Particle("C-14")], {"mass_numb": 13}, "", InvalidParticleError),
    ([Particle("Au 1+")], {"Z": 2}, "", InvalidParticleError),
    ([[]], {}, "", TypeError),
    (["D"], {}, ".recombine()", ChargeError),
    (["Fe 26+"], {}, ".ionize()", InvalidIonError),
    (["Fe 6+"], {}, ".ionize(-1)", ValueError),
    (["Fe 25+"], {}, ".recombine(0)", ValueError),
    (["Fe 6+"], {}, ".ionize(4.6)", TypeError),
    (["Fe 25+"], {}, ".recombine(8.2)", TypeError),
    (["e-"], {}, ".ionize()", InvalidElementError),
    (["e+"], {}, ".recombine()", InvalidElementError),
    (["e+"], {}, ".is_category('invalid_category')", ParticleError),
    (["e+"], {}, ".is_category(require='element', exclude='element')", ParticleError),
    (["H", 1], {}, "", TypeError),
    (["H", 1, 1], {}, "", TypeError),
    (["e-"], {}, ".nucleus", InvalidElementError),
]


@pytest.mark.parametrize(
    ("args", "kwargs", "attribute", "exception"), test_Particle_error_table
)
def test_Particle_errors(args, kwargs, attribute, exception) -> None:
    """
    Test that the appropriate exceptions are raised during the creation
    and use of a `~plasmapy.particles.Particle` object.
    """
    with pytest.raises(exception):
        exec(f"Particle(*args, **kwargs){attribute}")  # noqa: S102
        pytest.fail(
            f"The following command: "
            f"\n\n  {call_string(Particle, args, kwargs)}{attribute}\n\n"
            f"did not raise a {exception.__name__} as expected"
        )


# arg, kwargs, attribute, exception
test_Particle_warning_table = [
    ("H----", {}, "", ParticleWarning),
    ("alpha", {"mass_numb": 4}, "", ParticleWarning),
    ("alpha", {"Z": 2}, "", ParticleWarning),
]


@pytest.mark.parametrize(
    ("arg", "kwargs", "attribute", "warning"), test_Particle_warning_table
)
def test_Particle_warnings(arg, kwargs, attribute, warning) -> None:
    """
    Test that the appropriate warnings are issued during the creation
    and use of a `~plasmapy.particles.Particle` object.
    """
    with pytest.warns(warning) as record:
        exec(f"Particle(arg, **kwargs){attribute}")  # noqa: S102
        if not record:
            pytest.fail(
                f"The following command: "
                f"\n\n >>> {call_string(Particle, arg, kwargs)}{attribute}\n\n"
                f"did not issue a {warning.__name__} as expected"
            )


def test_Particle_cmp() -> None:
    """Test ``__eq__`` and ``__ne__`` in the Particle class."""
    proton1 = Particle("p+")
    proton2 = Particle("proton")
    electron = Particle("e-")

    assert proton1 == proton2, "Particle('p+') == Particle('proton') is False."
    assert proton1 != electron, "Particle('p+') == Particle('e-') is True."

    assert electron != 1

    assert electron != "dfasdf"


@pytest.mark.parametrize("particle", ["p+", "D+", "T+", "alpha"])
def test_particle_equality_special_nuclides(particle) -> None:
    particle_from_string = Particle(particle)
    particle_from_numbers = Particle(
        particle_from_string.element_name,
        Z=particle_from_string.charge_number,
        mass_numb=particle_from_string.mass_number,
    )
    assert particle_from_string == particle_from_numbers
    assert particle_from_string._attributes == particle_from_numbers._attributes


nuclide_mass_and_mass_equiv_table = [
    ("n", "neutron"),
    ("p+", "proton"),
    ("H-1", "p+"),
    ("H-1 0+", "p+"),
    ("D", "D+"),
    ("T", "T+"),
    ("He-4", "alpha"),
    ("Fe-56", "Fe-56 26+"),
]


@pytest.mark.parametrize(("isotope", "ion"), nuclide_mass_and_mass_equiv_table)
def test_particle_class_mass_nuclide_mass(isotope: str, ion: str) -> None:
    """
    Test that the ``mass`` and ``nuclide_mass`` attributes return
    equivalent values when appropriate.  The inputs should generally be
    an isotope with no charge information, and a fully ionized ion of
    that isotope, in order to make sure that the nuclide mass of the
    isotope equals the mass of the fully ionized ion.  This method may
    also check neutrons and protons.
    """

    Isotope = Particle(isotope)
    Ion = Particle(ion)

    if Isotope.categories & {"isotope", "baryon"} and Ion.categories & {
        "ion",
        "baryon",
    }:
        particle = Isotope.symbol

        assert Isotope.nuclide_mass == Ion.mass, (
            f"Particle({particle!r}).nuclide_mass does not equal "
            f"Particle({particle!r}).mass"
        )

    else:
        inputerrmsg = (
            f"isotope = {isotope!r} and ion = {ion!r} are "
            f"not valid inputs to this test. The inputs should be "
            f"an isotope with no charge information, and a fully "
            f"ionized ion of that isotope, in order to make sure "
            f"that the nuclide mass of the isotope equals the mass "
            f"of the ion."
        )

        assert Isotope.isotope and not Isotope.ion, inputerrmsg  # noqa: PT018
        assert Isotope.isotope == Ion.isotope, inputerrmsg
        assert Ion.charge_number == Ion.atomic_number, inputerrmsg

        assert Isotope.nuclide_mass == Ion.mass, (
            f"The nuclide mass of {isotope} does not equal the mass of {ion} "
            f"which is the fully ionized ion of that isotope. The results of "
            f"the test are:\n\n"
            f"Particle({ion!r}).mass = {Ion.mass}\n"
            f"Particle({isotope!r}).nuclide_mass = {Isotope.nuclide_mass}"
            "\n"
        )


@pytest.mark.slow
def test_particle_half_life_string() -> None:
    """
    Find the first isotope where the half-life is stored as a string
    (because the uncertainties are too great), and tests that requesting
    the half-life of that isotope causes a `MissingParticleDataWarning`
    whilst returning a string.
    """

    for isotope in known_isotopes():
        half_life = data_about_isotopes[isotope.isotope].get("half-life", None)
        if isinstance(half_life, str):
            break

    with pytest.warns(MissingParticleDataWarning):
        assert isinstance(Particle(isotope).half_life, str)


@pytest.mark.parametrize(
    ("p", "is_one"), [(Particle("e-"), True), (Particle("p+"), False)]
)
def test_particle_is_electron(p, is_one: bool) -> None:
    assert p.is_electron == is_one


def test_particle_bool_error() -> None:
    with pytest.raises(ParticleError):
        bool(Particle("e-"))


def test_particle_inversion(particle_antiparticle_pair) -> None:
    """Test that particles have the correct antiparticles."""
    particle, antiparticle = particle_antiparticle_pair
    assert particle.antiparticle == antiparticle, (
        f"The antiparticle of {particle} is found to be "
        f"{~Particle(particle)} instead of {antiparticle}."
    )


def test_antiparticle_inversion(particle_antiparticle_pair) -> None:
    """Test that antiparticles have the correct antiparticles."""
    particle, antiparticle = particle_antiparticle_pair
    assert antiparticle.antiparticle == particle, (
        f"The antiparticle of {antiparticle} is found to be "
        f"{~Particle(antiparticle)} instead of {particle}."
    )


def test_unary_operator_for_elements() -> None:
    with pytest.raises(ParticleError):
        Particle("C").antiparticle  # noqa: B018


class Test_antiparticle_properties_inversion:
    """
    Test particle and antiparticle inversion and properties for Particle
    instances.
    """

    def test_inverted_inversion(self, particle) -> None:
        """
        Test that the antiparticle of the antiparticle of a particle is
        the original particle.
        """
        assert particle == ~~particle, (
            f"~~{particle!r} equals {~~particle!r} instead of {particle!r}."
        )

    def test_opposite_charge(self, particle, opposite) -> None:
        """
        Test that a particle and its antiparticle have the opposite
        charge.
        """
        assert particle.charge_number == -opposite.charge_number, (
            f"The charges of {particle} and {opposite} are not "
            f"opposites, as expected of a particle/antiparticle pair."
        )

    def test_equal_mass(self, particle, opposite) -> None:
        """
        Test that a particle and its antiparticle have the same mass.
        """
        assert particle._attributes["mass"] == opposite._attributes["mass"], (
            f"The masses of {particle} and {opposite} are not equal, "
            f"as expected of a particle/antiparticle pair."
        )

    def test_antiparticle_attribute_and_operator(self, particle) -> None:
        """
        Test that the Particle.antiparticle attribute returns the same
        value as the unary ~ (invert) operator acting on the same
        Particle instance.
        """
        assert particle.antiparticle == ~particle, (
            f"{particle!r}.antiparticle returned "
            f"{particle.antiparticle}, whereas ~{particle!r} "
            f"returned {~particle}."
        )


@pytest.mark.parametrize("arg", ["e-", "D+", "Fe 25+", "H-", "mu+"])
def test_particleing_a_particle(arg) -> None:
    """
    Test that Particle(arg) is equal to Particle(Particle(arg)), but is
    not the same object in memory.
    """
    particle = Particle(arg)

    assert particle == Particle(particle), (
        f"Particle({arg!r}) does not equal Particle(Particle({arg!r})."
    )

    assert particle == Particle(Particle(Particle(particle))), (
        f"Particle({arg!r}) does not equal Particle(Particle(Particle({arg!r}))."
    )

    assert particle is not Particle(particle), (
        f"Particle({arg!r}) is the same object in memory as "
        f"Particle(Particle({arg!r})), when it is intended to "
        f"create a new object in memory (e.g., a copy)."
    )


@pytest.mark.parametrize(
    "key",
    [Particle("H"), Particle("e+"), CustomParticle(2 * 126.90447 * u.u, 0 * u.C, "I2")],
)
def test_that_object_can_be_dict_key(key):
    """
    Test that ``key`` can be used as the key of a `dict`.

    This test will fail if ``key`` does not equal itself or if ``key``
    is not hashable.  If this test fails, there is a problem with the
    ``__eq__`` and ``__hash__`` methods of ``key``.

    In most cases, objects that are mutable should not be hashable since
    they may change.

    """
    # TODO: I wrote this to be pretty general since I felt like
    # TODO: procrastinating other things, so we can probably put this
    # TODO: into utils.pytest_helpers later on.  There are likely other
    # TODO: classes that should be able to be used as keys of dicts.

    value = 42

    try:
        dictionary = {key: value}
    except Exception as exc:
        error_message = f"{key} is not a valid key for a dict. "
        if not isinstance(key, collections.abc.Hashable):
            error_message += f"{key} is not hashable. "
        try:
            key_equals_itself = key == key  # noqa: PLR0124
        except Exception:  # noqa: BLE001
            error_message += f"{key} == {key} cannot be evaluated. "
        else:
            if not key_equals_itself:
                error_message += f"{key} does not equal itself."
        raise TypeError(error_message) from exc

    assert dictionary[key] == value


customized_particle_tests = [
    (DimensionlessParticle, {"mass": 1.0, "charge": -1.0}, "mass", 1.0),
    (DimensionlessParticle, {"mass": 0.0, "charge": -1.0}, "charge", -1.0),
    (DimensionlessParticle, {}, "mass", np.nan),
    (DimensionlessParticle, {}, "charge", np.nan),
    (DimensionlessParticle, {"mass": np.inf}, "mass", np.inf),
    (DimensionlessParticle, {"charge": np.inf}, "charge", np.inf),
    (DimensionlessParticle, {"charge": 1.0 * u.dimensionless_unscaled}, "charge", 1.0),
    (DimensionlessParticle, {"mass": 1.0 * u.dimensionless_unscaled}, "mass", 1.0),
    (CustomParticle, {}, "mass", np.nan * u.kg),
    (CustomParticle, {}, "charge", np.nan * u.C),
    (CustomParticle, {"mass": 1.1 * u.kg, "charge": -0.1 * u.C}, "mass", 1.1 * u.kg),
    (CustomParticle, {"charge": -0.1 * u.C}, "charge", -0.1 * u.C),
    (CustomParticle, {"mass": np.inf * u.g}, "mass", np.inf * u.kg),
    (CustomParticle, {"mass": "100.0 g"}, "mass", 100.0 * u.g),
    (CustomParticle, {"charge": -np.inf * u.kC}, "charge", -np.inf * u.C),
    (CustomParticle, {"charge": "5.0 C"}, "charge", 5.0 * u.C),
    (CustomParticle, {"Z": 1}, "charge", const.e.si),
    (CustomParticle, {"Z": 1.5}, "charge", 1.5 * const.e.si),
    (CustomParticle, {"Z": 1.5}, "charge_number", 1.5),
    (CustomParticle, {"charge": 1.3 * const.e.si}, "charge_number", 1.3),
]


@pytest.mark.parametrize(
    ("cls", "kwargs", "attr", "expected"), customized_particle_tests
)
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_custom_particles(cls, kwargs, attr, expected) -> None:
    """Test the attributes of dimensionless and custom particles."""
    instance = cls(**kwargs)
    value = getattr(instance, attr)
    if not u.isclose(value, expected, equal_nan=True):
        pytest.fail(
            f"{call_string(cls, kwargs=kwargs)}.{attr} should return a value "
            f"of {expected}, but instead returned a value of {value}."
        )


@pytest.mark.parametrize(
    ("cls", "symbol", "expected"),
    [
        (CustomParticle, None, "CustomParticle(mass=nan kg, charge=nan C)"),
        (CustomParticle, "η", "η"),
        (DimensionlessParticle, None, "DimensionlessParticle(mass=nan, charge=nan)"),
        (DimensionlessParticle, "η", "η"),
    ],
)
def test_custom_particle_symbol(cls, symbol, expected) -> None:
    instance = cls(symbol=symbol)
    assert instance.symbol == expected


custom_particle_categories_table = [
    ({"charge": 0.0 * u.C}, {"custom", "uncharged"}),
    ({"charge": 1.0 * u.C}, {"custom", "charged"}),
    ({}, {"custom"}),
]


@pytest.mark.parametrize(("kwargs", "expected"), custom_particle_categories_table)
def test_custom_particle_categories(kwargs, expected) -> None:
    """Test that CustomParticle.categories behaves as expected."""
    custom_particle = CustomParticle(**kwargs)
    assert custom_particle.categories == expected


custom_particle_is_category_table = [
    ({"charge": 0 * u.C}, {"require": "charged"}, False),
    ({"charge": 0 * u.C}, {"exclude": "charged"}, True),
    ({"charge": 0 * u.C}, {"require": "uncharged"}, True),
    ({"charge": 0 * u.C}, {"exclude": "uncharged"}, False),
    ({"charge": 1 * u.C}, {"require": "charged"}, True),
    ({"charge": 1 * u.C}, {"exclude": "charged"}, False),
    ({"charge": 1 * u.C}, {"require": "uncharged"}, False),
    ({"charge": 1 * u.C}, {"exclude": "uncharged"}, True),
    ({}, {"any_of": {"charged", "uncharged"}}, False),
]


@pytest.mark.parametrize(
    ("kwargs_to_custom_particle", "kwargs_to_is_category", "expected"),
    custom_particle_is_category_table,
)
def test_custom_particle_is_category(
    kwargs_to_custom_particle,
    kwargs_to_is_category,
    expected,
) -> None:
    """Test that CustomParticle.is_category works as expected."""
    custom_particle = CustomParticle(**kwargs_to_custom_particle)
    actual = custom_particle.is_category(**kwargs_to_is_category)
    assert actual == expected


custom_particle_errors = [
    (DimensionlessParticle, {"mass": -1e-36}, InvalidParticleError),
    (DimensionlessParticle, {"mass": [1, 1]}, InvalidParticleError),
    (DimensionlessParticle, {"charge": [-1, 1]}, InvalidParticleError),
    (DimensionlessParticle, {"mass": True}, InvalidParticleError),
    (
        DimensionlessParticle,
        {"mass": np.array([1, 2]) * u.dimensionless_unscaled},
        InvalidParticleError,
    ),
    (
        DimensionlessParticle,
        {"charge": np.array([1, 2]) * u.dimensionless_unscaled},
        InvalidParticleError,
    ),
    (
        DimensionlessParticle,
        {"charge": (1 + 2j) * u.dimensionless_unscaled},
        InvalidParticleError,
    ),
    (CustomParticle, {"charge": 5 + 2j}, InvalidParticleError),
    (CustomParticle, {"mass": 5 + 2j}, InvalidParticleError),
    (CustomParticle, {"charge": np.complex128(5 + 2j)}, InvalidParticleError),
    (CustomParticle, {"mass": np.complex128(5 + 2j)}, InvalidParticleError),
    (CustomParticle, {"mass": -1e-36 * u.kg}, InvalidParticleError),
    (CustomParticle, {"mass": "not a mass"}, InvalidParticleError),
    (CustomParticle, {"mass": "5.0 km"}, InvalidParticleError),
    (CustomParticle, {"mass": np.array([1, 1]) * u.kg}, InvalidParticleError),
    (CustomParticle, {"charge": np.array([1, 1]) * u.C}, InvalidParticleError),
    (CustomParticle, {"charge": (5 + 2j) * u.C}, InvalidParticleError),
    (CustomParticle, {"mass": (5 + 2j) * u.kg}, InvalidParticleError),
    (CustomParticle, {"charge": np.complex128(5 + 2j) * u.C}, InvalidParticleError),
    (CustomParticle, {"mass": np.complex128(5 + 2j) * u.kg}, InvalidParticleError),
    (CustomParticle, {"charge": "not a charge"}, InvalidParticleError),
    (CustomParticle, {"charge": "5.0 km"}, InvalidParticleError),
    (CustomParticle, {"charge": 1 * u.C, "Z": -1}, InvalidParticleError),
    (CustomParticle, {"charge": 1}, InvalidParticleError),
]


@pytest.mark.parametrize(("cls", "kwargs", "exception"), custom_particle_errors)
def test_customized_particles_errors(cls, kwargs, exception) -> None:
    """
    Test that attempting to create invalid dimensionless or custom particles
    results in an InvalidParticleError.
    """
    with pytest.raises(exception):
        cls(**kwargs)


customized_particle_repr_table = [
    (
        CustomParticle,
        {"mass": 5.12 * u.kg, "charge": 6.2 * u.C},
        "CustomParticle(mass=5.12 kg, charge=6.2 C)",
    ),
    (
        DimensionlessParticle,
        {"mass": 5.2, "charge": 6.3},
        "DimensionlessParticle(mass=5.2, charge=6.3)",
    ),
]


@pytest.mark.parametrize(
    ("cls", "kwargs", "expected_repr"), customized_particle_repr_table
)
def test_customized_particle_repr(cls, kwargs, expected_repr) -> None:
    """Test the string representations of dimensionless and custom particles."""
    instance = cls(**kwargs)
    from_repr = repr(instance)
    from_str = str(instance)
    if not expected_repr == from_repr == from_str:
        pytest.fail(
            f"Problem with a string representation of {cls.__name__} "
            f"with kwargs = {kwargs}.\n\n"
            f"expected_repr = {expected_repr}"
            f"from_str: {from_str}"
            f"from_repr: {from_repr}"
        )


custom_particle_str_table = [
    (
        {"mass": 5.12 * u.kg, "charge": 6.2 * u.C},
        "CustomParticle(mass=5.12 kg, charge=6.2 C)",
    ),
    (
        {"mass": 5.12 * u.kg, "charge": 6.2 * u.C, "symbol": "I2"},
        "I2",
    ),
]


@pytest.mark.parametrize(("kwargs", "expected_str"), custom_particle_str_table)
def test_custom_particle_str(kwargs, expected_str) -> None:
    """Test the string representation of a custom particle."""
    instance = CustomParticle(**kwargs)
    assert str(instance) == expected_str


@pytest.mark.parametrize("cls", [CustomParticle, DimensionlessParticle])
@pytest.mark.parametrize("not_a_str", [1, u.kg])
def test_typeerror_redefining_symbol(cls, not_a_str) -> None:
    """Test that the symbol attribute cannot be set to something besides a string"""
    instance = cls()
    with pytest.raises(TypeError):
        instance.symbol = not_a_str


custom_particles_from_json_tests = [
    (
        DimensionlessParticle,
        {"mass": 5.2, "charge": 6.3, "symbol": "ξ"},
        '{"plasmapy_particle": {"type": "DimensionlessParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": { \
            "args": [], \
            "kwargs": {"mass": 5.2, "charge": 6.3, "symbol": "ξ"}}}}',
        None,
    ),
    (
        DimensionlessParticle,
        {"mass": 5.2, "symbol": "ξ"},
        '{"plasmapy_particle": {"type": "DimensionlessParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": { \
            "args": [], \
            "kwargs": {"mass": 5.2, "charge": NaN, "symbol": "ξ"}}}}',
        None,
    ),
    (
        CustomParticle,
        {"mass": 5.12 * u.kg, "charge": 6.2 * u.C, "symbol": "ξ"},
        '{"plasmapy_particle": {"type": "CustomParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": { \
            "args": [], \
            "kwargs": {"mass": "5.12 kg", "charge": "6.2 C", '
        '"symbol": "ξ"}}}}',
        None,
    ),
    (
        CustomParticle,
        {"mass": 5.12 * u.kg},
        '{"plasmapy_particle": {"type": "CustomParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": { \
            "args": [], \
            "kwargs": {"mass": "5.12 kg", "charge": "nan C"}}}}',
        None,
    ),
    (
        DimensionlessParticle,
        {"mass": 5.2, "charge": 6.3},
        '{"plasmapy_particle": {"notatype": "DimensionlessParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": { \
            "args": [], \
            "kwargs": {"mass": 5.2, "charge": 6.3}}}}',
        InvalidElementError,
    ),
    (
        CustomParticle,
        {"mass": 5.12 * u.kg},
        '{"plasmapy_particle": {"notatype": "CustomParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": { \
            "args": [], \
            "kwargs": {"mass": "5.12 kg", "charge": "nan C"}}}}',
        InvalidElementError,
    ),
    (
        DimensionlessParticle,
        {"mass": 5.2, "charge": 6.3},
        '{"plasmapy_particle": {"type": "DimensionlessParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "fake__init__": { \
            "args": [], \
            "kwargs": {"mass": 5.2, "charge": 6.3}}}}',
        InvalidElementError,
    ),
    (
        CustomParticle,
        {"mass": 5.12 * u.kg},
        '{"plasmapy_particle": {"type": "CustomParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "fake__init__": { \
            "args": [], \
            "kwargs": {"mass": "5.12 kg", "charge": "nan C"}}}}',
        InvalidElementError,
    ),
]


@pytest.mark.parametrize(
    ("cls", "kwargs", "json_string", "expected_exception"),
    custom_particles_from_json_tests,
)
def test_custom_particles_from_json_string(
    cls, kwargs, json_string: str, expected_exception
) -> None:
    """Test the attributes of dimensionless and custom particles generated from
    JSON representation"""
    if expected_exception is None:
        instance = cls(**kwargs)
        instance_from_json = json_loads_particle(json_string)
        assert u.isclose(instance.mass, instance_from_json.mass, equal_nan=True), (
            pytest.fail(
                f"Expected a mass value of {instance.mass}\n"
                f"Received a mass value of {instance_from_json.mass}"
            )
        )
        assert u.isclose(instance.charge, instance_from_json.charge, equal_nan=True), (
            pytest.fail(
                f"Expected a charge value of {instance.charge}\n"
                f"Received a charge value of {instance_from_json.charge}"
            )
        )
    else:
        with pytest.raises(expected_exception):
            instance_from_json = json_loads_particle(json_string)
            pytest.fail(
                f"{cls.__name__} with ({json_string})"
                f" did not raise: {expected_exception.__name__}."
            )


@pytest.mark.parametrize(
    ("cls", "kwargs", "json_string", "expected_exception"),
    custom_particles_from_json_tests,
)
def test_custom_particles_from_json_file(
    cls, kwargs, json_string: str, expected_exception
) -> None:
    """Test the attributes of dimensionless and custom particles generated from
    JSON representation"""
    if expected_exception is None:
        instance = cls(**kwargs)
        test_file_object = io.StringIO(json_string)
        instance_from_json = json_load_particle(test_file_object)
        assert u.isclose(instance.mass, instance_from_json.mass, equal_nan=True), (
            pytest.fail(
                f"Expected a mass value of {instance.mass}\n"
                f"Received a mass value of {instance_from_json.mass}"
            )
        )
        assert u.isclose(instance.charge, instance_from_json.charge, equal_nan=True), (
            pytest.fail(
                f"Expected a charge value of {instance.charge}\n"
                f"Received a charge value of {instance_from_json.charge}"
            )
        )
    else:
        with pytest.raises(expected_exception):
            test_file_object = io.StringIO(json_string)
            instance_from_json = json_load_particle(test_file_object)
            pytest.fail(
                f"{cls.__name__} with ({json_string})"
                f" did not raise: {expected_exception.__name__}."
            )


particles_from_json_tests = [
    (
        Particle,
        {"argument": "Pb"},
        '{"plasmapy_particle": {"type": "Particle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": {"args": ["Pb"], "kwargs": {}}}}',
        None,
    ),
    (
        Particle,
        {"argument": "e-"},
        '{"plasmapy_particle": {"type": "Particle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": {"args": ["e-"], "kwargs": {}}}}',
        None,
    ),
    (
        Particle,
        {"argument": "e-"},
        '{"plasmapy_particle": {"notatype": "Particle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": {"args": ["e-"], "kwargs": {}}}}',
        InvalidElementError,
    ),
    (
        Particle,
        {"argument": "e-"},
        '{"plasmapy_particle": {"type": "Particle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "fake__init__": {"args": ["e-"], "kwargs": {}}}}',
        InvalidElementError,
    ),
]


@pytest.mark.parametrize(
    ("cls", "kwargs", "json_string", "expected_exception"), particles_from_json_tests
)
def test_particles_from_json_string(
    cls, kwargs, json_string: str, expected_exception
) -> None:
    """Test the attributes of Particle objects created from JSON representation."""
    if expected_exception is None:
        instance = cls(**kwargs)
        instance_from_json = json_loads_particle(json_string)
        expected_particle = instance.symbol
        actual_particle = instance_from_json.symbol
        assert expected_particle == actual_particle, pytest.fail(
            f"Expected {expected_particle}\nGot {actual_particle}"
        )
    else:
        with pytest.raises(expected_exception):
            instance_from_json = json_loads_particle(json_string)
            pytest.fail(
                f"{cls.__name__} with ({json_string})"
                f" did not raise: {expected_exception.__name__}."
            )


@pytest.mark.parametrize(
    ("cls", "kwargs", "json_string", "expected_exception"), particles_from_json_tests
)
def test_particles_from_json_file(
    cls, kwargs, json_string: str, expected_exception
) -> None:
    """Test the attributes of Particle objects created from JSON representation."""
    if expected_exception is None:
        instance = cls(**kwargs)
        test_file_object = io.StringIO(json_string)
        test_file_object.seek(0, io.SEEK_SET)
        instance_from_json = json_load_particle(test_file_object)
        expected_particle = instance.symbol
        actual_particle = instance_from_json.symbol
        assert expected_particle == actual_particle, pytest.fail(
            f"Expected {expected_particle}\nGot {actual_particle}"
        )
    else:
        with pytest.raises(expected_exception):
            test_file_object = io.StringIO(json_string)
            instance_from_json = json_load_particle(test_file_object)
            pytest.fail(
                f"{cls.__name__} with ({json_string})"
                f" did not raise: {expected_exception.__name__}."
            )


particle_json_repr_table = [
    (
        Particle,
        {"argument": "lead"},
        '{"plasmapy_particle": {"type": "Particle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", \
        "__init__": {"args": ["Pb"], "kwargs": {}}}}',
    ),
    (
        Particle,
        {"argument": "lead"},
        '{"plasmapy_particle": {"type": "Particle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": {"args": ["Pb"], "kwargs": {}}}}',
    ),
    (
        CustomParticle,
        {"mass": 5.12 * u.kg, "charge": 6.2e-19 * u.C, "symbol": "ξ"},
        '{"plasmapy_particle": {"type": "CustomParticle", \
        "module": "plasmapy.particles.particle_class", \
        "date_created": "...", "__init__": {\
            "args": [], \
            "kwargs": {"mass": "5.12 kg", "charge": "6.2e-19 C", \
            "charge_number": "3.869735626165673", "symbol": "ξ"}}}}',
    ),
    (
        DimensionlessParticle,
        {"mass": 5.2, "charge": 6.3, "symbol": "ξ"},
        '{"plasmapy_particle": {"type": "DimensionlessParticle",\
        "module": "plasmapy.particles.particle_class",\
        "date_created": "...", "__init__": {\
            "args": [], \
            "kwargs": {"mass": 5.2, "charge": 6.3, "symbol": "ξ"}}}}',
    ),
]


@pytest.mark.parametrize(("cls", "kwargs", "expected_repr"), particle_json_repr_table)
def test_particle_to_json_string(cls, kwargs, expected_repr) -> None:
    """Test the JSON representations of normal, dimensionless and custom particles."""
    instance = cls(**kwargs)
    json_repr = instance.json_dumps()
    test_dict = json.loads(json_repr)["plasmapy_particle"]
    expected_repr = json.loads(expected_repr)["plasmapy_particle"]
    assert test_dict["type"] == expected_repr["type"], pytest.fail(
        f"Problem with JSON representation of {cls.__name__} "
        f"with kwargs = {kwargs}.\n\n"
        f"expected type = {expected_repr['type']}\n\n"
        f"got type: {test_dict['type']}"
    )
    assert expected_repr["__init__"] == test_dict["__init__"], pytest.fail(
        f"Problem with JSON representation of {cls.__name__} "
        f"with kwargs = {kwargs}.\n\n"
        f"expected_repr = {expected_repr['__init__']}.\n\n"
        f"json_repr: {test_dict['__init__']}"
    )


@pytest.mark.parametrize(("cls", "kwargs", "expected_repr"), particle_json_repr_table)
def test_particle_to_json_file(cls, kwargs, expected_repr) -> None:
    """Test the JSON representations of normal, dimensionless and custom particles."""
    instance = cls(**kwargs)
    test_file_object = io.StringIO("")
    instance.json_dump(test_file_object)
    test_file_object.seek(0, io.SEEK_SET)
    json_repr = test_file_object.read()
    test_dict = json.loads(json_repr)["plasmapy_particle"]
    expected_repr = json.loads(expected_repr)["plasmapy_particle"]
    assert test_dict["type"] == expected_repr["type"], pytest.fail(
        f"Problem with JSON representation of {cls.__name__} "
        f"with kwargs = {kwargs}.\n\n"
        f"expected type = {expected_repr['type']}\n\n"
        f"got type: {test_dict['type']}"
    )
    assert expected_repr["__init__"] == test_dict["__init__"], pytest.fail(
        f"Problem with JSON representation of {cls.__name__} "
        f"with kwargs = {kwargs}.\n\n"
        f"expected_repr = {expected_repr['__init__']}.\n\n"
        f"json_repr: {test_dict['__init__']}"
    )


def test_particle_is_category_valid_categories() -> None:
    """Test the location where valid categories may be accessed."""
    some_valid_categories = {
        "charged",
        "custom",
        "electron",
        "fermion",
        "ion",
        "isotope",
        "lepton",
        "matter",
        "nonmetal",
        "uncharged",
    }
    assert some_valid_categories.issubset(valid_categories)


def test_CustomParticle_cmp() -> None:
    """Test ``__eq__`` and ``__ne__`` in the CustomParticle class."""
    particle1 = CustomParticle(2 * 126.90447 * u.u, 0 * u.C, "I2")
    particle2 = CustomParticle(2 * 126.90447 * u.u, 0 * u.C, "I2")
    other = CustomParticle(2 * 126.90447 * u.u, e.si, "I2 +")

    assert particle1 == particle2, (
        "CustomParticle instances that should be equal are not."
    )
    assert particle1 != other, "CustomParticle instances should not be equal, but are."

    assert particle1 != 1


@pytest.mark.parametrize(
    ("attr", "arg"),
    [
        ("charge", 2 * u.C),
        ("mass", 2 * u.kg),
        ("charge_number", 2),
        ("symbol", "😺"),
    ],
)
def test_CustomParticle_setters(attr, arg) -> None:
    custom_particle = CustomParticle()
    setattr(custom_particle, attr, arg)
    output = getattr(custom_particle, attr)
    assert arg == output


test_molecule_table = [
    (2 * 126.90447 * u.u, 0 * u.C, "I2", "I2", None),
    (2 * 126.90447 * u.u, e.si, "I2 1+", "I2 1+", None),
    (2 * 126.90447 * u.u, e.si, "I2 1+", "I2", 1),
    (2 * 126.90447 * u.u, e.si, "II 1+", "II", 1),
]


@pytest.mark.parametrize(
    ("args", "kwargs", "expected"),
    [
        ([], {}, CustomParticle()),
        ([1 * u.kg], {}, CustomParticle(mass=1 * u.kg)),
        ([2 * u.C], {}, CustomParticle(charge=2 * u.C)),
        ([3 * u.kg, 4 * u.C], {}, CustomParticle(mass=3 * u.kg, charge=4 * u.C)),
        ([5 * u.C, 6 * u.kg], {}, CustomParticle(mass=6 * u.kg, charge=5 * u.C)),
        ([7 * u.kg], {"Z": 8.9}, CustomParticle(mass=7 * u.kg, Z=8.9)),
        ([], {"symbol": "..."}, CustomParticle(symbol="...")),
    ],
)
def test_CustomParticle_from_quantities(args, kwargs, expected) -> None:
    actual = CustomParticle._from_quantities(*args, **kwargs)
    assert actual == expected


@pytest.mark.parametrize(
    ("args", "kwargs", "exception"),
    [
        ([1 * u.C], {"Z": 2}, InvalidParticleError),
        ([3 * u.m], {}, InvalidParticleError),
        ([4 * u.kg, "invalid"], {}, InvalidParticleError),
    ],
)
def test_CustomParticle_from_quantities_errors(args, kwargs, exception) -> None:
    with pytest.raises(exception):
        CustomParticle._from_quantities(*args, **kwargs)


@pytest.mark.parametrize(("m", "Z", "symbol", "m_symbol", "m_Z"), test_molecule_table)
def test_molecule(m, Z, symbol, m_symbol, m_Z) -> None:
    """Test ``molecule`` function."""
    assert CustomParticle(m, Z, symbol) == molecule(m_symbol, m_Z)


test_molecule_error_table = [
    ("Zz", None),
    ("", None),
    ("I2+", 2),
    ("Iii", None),
    ("e2H3", None),
]


@pytest.mark.parametrize(("symbol", "Z"), test_molecule_error_table)
def test_molecule_error(symbol, Z) -> None:
    """Test the error raised in case of a bad molecule symbol."""
    with pytest.raises(InvalidParticleError):
        molecule(symbol, Z)


def test_molecule_other() -> None:
    """Test fallback to |Particle| object and warning in case of redundant charge."""
    assert Particle("I") == molecule("I")

    with pytest.warns(ParticleWarning):
        assert CustomParticle(2 * 126.90447 * u.u, e.si, "I2 1+") == molecule(
            "I2 1+", Z=1
        )


def test_undefined_charge() -> None:
    """Test that a particle with an undefined charge returns |nan| C."""
    H_particle = Particle("H")
    assert u.isclose(H_particle.charge, np.nan * u.C, equal_nan=True)


def test_undefined_standard_atomic_weight() -> None:
    """Test that a particle with an undefined standard atomic weight returns |nan| kg."""
    Pm_particle = Particle("Pm")
    assert u.isclose(Pm_particle.standard_atomic_weight, np.nan * u.kg, equal_nan=True)


def test_undefined_mass() -> None:
    """Test that a particle with an undefined mass returns |nan| kg."""
    tau_neutrino_particle = Particle("tau neutrino")
    assert u.isclose(tau_neutrino_particle.mass, np.nan * u.kg, equal_nan=True)


def test_undefined_mass_energy() -> None:
    """Test that a particle with an undefined mass energy returns |nan| J."""
    nu_tau_particle = Particle("nu_tau")
    assert u.isclose(nu_tau_particle.mass_energy, np.nan * u.J, equal_nan=True)


def test_isotopes_have_equivalent_ionization_energy() -> None:
    """Test that isotopes of the same element with the same charge have the same ionization energy."""
    U = Particle("U +1")
    U_235 = Particle("U-235 +1")
    U_238 = Particle("U-238 +1")

    assert U.ionization_energy == U_235.ionization_energy == U_238.ionization_energy


def test_different_ions_have_different_ionization_energy() -> None:
    """Test that isotopes of the same element with different charges have different ionization energies."""
    U_1 = Particle("U +1")
    U_2 = Particle("U +2")
    U_3 = Particle("U +3")

    assert U_1.ionization_energy != U_2.ionization_energy
    assert U_1.ionization_energy != U_3.ionization_energy
    assert U_2.ionization_energy != U_3.ionization_energy


def test_zero_charge_ionization_energy() -> None:
    """Test that the ionization energy is the same if initialized with a zero charge or just the element symbol."""
    C_0 = Particle("C +0")
    C = Particle("C")
    assert C_0.ionization_energy == C.ionization_energy


def test_deuterium_ionization_energy() -> None:
    """Test that the ionization energy of deuterium is the same as the ionization energy of H-2."""
    D = Particle("D")
    H_2 = Particle("H-2")
    H = Particle("H")
    assert D.ionization_energy == H_2.ionization_energy
    assert D.ionization_energy != H.ionization_energy


@pytest.mark.parametrize(
    ("particle_symbol", "expected_ionization_energy"),
    [
        ("H", (13.598434599702 * u.eV).to(u.J)),
        ("He", u.Quantity((24.587389011 * u.eV).to(u.J).value, u.J)),
        ("Li", u.Quantity((5.391714996 * u.eV).to(u.J).value, u.J)),
    ],
)
def test_particle_ionization_energy(
    particle_symbol, expected_ionization_energy
) -> None:
    particle = Particle(particle_symbol)
    assert u.isclose(
        particle.ionization_energy, expected_ionization_energy, rtol=1e-4
    ), f"Expected {expected_ionization_energy}, got {particle.ionization_energy}"


def test_undefined_ionization_energy() -> None:
    particle = Particle("tau neutrino")
    try:
        energy = particle.ionization_energy
        pytest.fail(f"Expected MissingParticleDataError, got {energy}")
    except MissingParticleDataError:
        pass


def test_undefined_electron_binding_energy() -> None:
    particle = Particle("tau neutrino")
    try:
        energy = particle.electron_binding_energy
        pytest.fail(f"Expected MissingParticleDataError, got {energy}")
    except MissingParticleDataError:
        pass


def test_warning_on_use_of_binding_energy() -> None:
    with pytest.warns(FutureWarning):
        particle = Particle("n")
        assert particle.binding_energy == particle.nuclear_binding_energy


def test_deuterium_electron_binding_energy() -> None:
    D = Particle("D")
    H_2 = Particle("H-2")
    H = Particle("H")
    assert D.electron_binding_energy == H_2.electron_binding_energy
    assert D.electron_binding_energy != H.electron_binding_energy


def test_other_isotopes_electron_binding_energy() -> None:
    C_12 = Particle("C-12")
    C_13 = Particle("C-13")
    C_14 = Particle("C-14")
    assert C_12.electron_binding_energy == C_13.electron_binding_energy
    assert C_12.electron_binding_energy == C_14.electron_binding_energy
    assert C_13.electron_binding_energy == C_14.electron_binding_energy


def test_isotope_ion_electron_binding_energy() -> None:
    C_12 = Particle("C-12 +1")
    C_13 = Particle("C-13 +1")
    C_14 = Particle("C-14 +1")
    assert C_12.electron_binding_energy == C_13.electron_binding_energy
    assert C_12.electron_binding_energy == C_14.electron_binding_energy
    assert C_13.electron_binding_energy == C_14.electron_binding_energy


def test_infinite_ionization() -> None:
    helium = Particle("He-4 0+")
    h_nucleus = helium.ionize(n=np.inf)
    helium.ionize(n=np.inf, inplace=True)
    assert h_nucleus == helium
    assert h_nucleus == Particle("He-4 2+")
    helium = Particle("He-4 0+")
    pytest.raises(TypeError, helium.ionize, n=0.5)
    pytest.raises(ValueError, helium.ionize, n=-1)
