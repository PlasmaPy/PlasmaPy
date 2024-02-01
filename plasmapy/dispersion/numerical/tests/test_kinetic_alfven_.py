"""Test functionality of Stix in `plasmapy.dispersion.numerical.kinetic_alfven_`."""
import astropy.units as u
import numpy as np
import pytest
from astropy.constants.si import c

from plasmapy.dispersion.numerical.kinetic_alfven_ import kinetic_alfven
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import PhysicsWarning

c_si_unitless = c.value


class TestKinetic_Alfven:
    _kwargs_single_valued = {
        "B": 8.3e-9 * u.T,
        "ion": Particle("p+"),
        "k": np.logspace(-7, -2, 2) * u.rad / u.m,
        "n_i": 5 * u.m**-3,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "theta": 30 * u.deg,
        "gamma_e": 3,
        "gamma_i": 3,
        "Z": 1,
    }

    @pytest.mark.parametrize(
        ("kwargs", "_error"),
        [
            ({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "B": [1e-9, 2e-9, 3e-9] * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": -1 * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": 5 * u.m}, u.UnitTypeError),
            ({**_kwargs_single_valued, "ion": "not a particle"}, InvalidParticleError),
            ({**_kwargs_single_valued, "k": np.ones((3, 2)) * u.rad / u.m}, ValueError),
            ({**_kwargs_single_valued, "k": 0 * u.rad / u.m}, ValueError),
            ({**_kwargs_single_valued, "k": -1.0 * u.rad / u.m}, ValueError),
            ({**_kwargs_single_valued, "k": 5 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "n_i": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "n_i": [5e6, 6e6] * u.m**-3}, ValueError),
            ({**_kwargs_single_valued, "n_i": -5e6 * u.m**-3}, ValueError),
            ({**_kwargs_single_valued, "n_i": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_e": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "T_e": [1.4e6, 1.7e6] * u.K}, ValueError),
            ({**_kwargs_single_valued, "T_e": -10 * u.eV}, ValueError),
            ({**_kwargs_single_valued, "T_e": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_i": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "T_i": [4e5, 5e5] * u.K}, ValueError),
            ({**_kwargs_single_valued, "T_i": -1 * u.eV}, ValueError),
            ({**_kwargs_single_valued, "T_i": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            ({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError),
            ({**_kwargs_single_valued, "gamma_e": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "gamma_i": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "Z": "wrong type"}, TypeError),
        ],
    )
    def test_raises(self, kwargs, _error) -> None:
        """Test scenarios that raise an `Exception`."""
        with pytest.raises(_error):
            kinetic_alfven(**kwargs)

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            ({**_kwargs_single_valued, "theta": 0 * u.deg}, {"shape": (2,)}),
            (
                {
                    **_kwargs_single_valued,
                    "theta": 0 * u.deg,
                    "k": [1, 2, 3] * u.rad / u.m,
                },
                {"shape": (3,)},
            ),
            (
                {
                    **_kwargs_single_valued,
                    "theta": [10, 20, 30, 40, 50] * u.deg,
                    "k": [1, 2, 3] * u.rad / u.m,
                },
                {"shape": (3,)},
            ),
            (
                {**_kwargs_single_valued, "theta": [10, 20, 30, 40, 50] * u.deg},
                {"shape": (2,)},
            ),
        ],
    )
    @pytest.mark.filterwarnings("ignore::plasmapy.utils.exceptions.PhysicsWarning")
    def test_return_structure(self, kwargs, expected) -> None:
        """Test the structure of the returned values."""
        ws = kinetic_alfven(**kwargs)

        assert isinstance(ws, dict)

        for val in ws.values():
            assert isinstance(val, u.Quantity)
            assert val.unit == u.rad / u.s
            assert val.shape == expected["shape"]

    @pytest.mark.parametrize(
        ("kwargs", "_warning"),
        [
            # w/vT min PhysicsWarning
            (
                {
                    "k": 10 * u.rad / u.m,
                    "theta": 88 * u.deg,
                    "n_i": 0.05 * u.cm**-3,
                    "B": 2.2e-8 * u.T,
                    "T_e": 1.6e6 * u.K,
                    "T_i": 4.0e5 * u.K,
                    "ion": Particle("p+"),
                },
                PhysicsWarning,
            ),
            # w/vT max PhysicsWarning
            (
                {
                    "k": 0.000001 * u.rad / u.m,
                    "theta": 88 * u.deg,
                    "n_i": 0.05 * u.cm**-3,
                    "B": 2.2e-8 * u.T,
                    "T_e": 1.6e6 * u.K,
                    "T_i": 4.0e5 * u.K,
                    "ion": Particle("p+"),
                },
                PhysicsWarning,
            ),
            # w << w_ci PhysicsWarning
            (
                {
                    "k": 10e-8 * u.rad / u.m,
                    "theta": 88 * u.deg,
                    "n_i": 5 * u.cm**-3,
                    "B": 6.98e-8 * u.T,
                    "T_e": 1.6e6 * u.K,
                    "T_i": 4.0e5 * u.K,
                    "ion": Particle("p+"),
                },
                PhysicsWarning,
            ),
        ],
    )
    def test_warning(self, kwargs, _warning) -> None:
        """Test scenarios that raise a `Warning`."""
        with pytest.warns(_warning):
            kinetic_alfven(**kwargs)
