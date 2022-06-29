"""Test functionality of Stix in `plasmapy.dispersion.analytical.stix_`."""
import numpy as np
import pytest

from astropy import units as u
from astropy.constants.si import c

from plasmapy.dispersion.analytical.stix_ import stix
from plasmapy.formulary import gyrofrequency, plasma_frequency
from plasmapy.particles import Particle, ParticleList
from plasmapy.particles.exceptions import InvalidParticleError

c_si_unitless = c.value


class TestStix:
    _kwargs_single_valued = {
        "B": 8.3e-9 * u.T,
        "w": 0.001 * u.rad / u.s,
        "ions": [Particle("He+"), Particle("H+")],
        "n_i": [4.0e5, 2.0e5] * u.m**-3,
        "theta": 30 * u.deg,
    }

    @staticmethod
    def spd(B, w, ions, n_i):
        w = w.to(u.rad / u.s).value
        n_i = n_i.to(u.m**-3).value

        species = ions + [Particle("e-")]
        densities = np.zeros(n_i.size + 1)
        densities[:-1] = n_i
        densities[-1] = np.sum(n_i * ions.charge_number)

        # Generate the plasma parameters needed
        wps = []
        wcs = []
        for par, dens in zip(species, densities):
            wps.append(plasma_frequency(n=dens * u.m**-3, particle=par).value)
            wcs.append(gyrofrequency(B=B, particle=par, signed=True).value)

        # Stix method implemented
        S = np.ones_like(w, dtype=np.float64)
        P = np.ones_like(S)
        D = np.zeros_like(S)
        for wc, wp in zip(wcs, wps):
            S -= (wp**2) / (w**2 - wc**2)
            P -= (wp / w) ** 2
            D += ((wp**2) / (w**2 - wc**2)) * (wc / w)

        return S, P, D

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            ({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": -1 * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": 5 * u.m}, u.UnitTypeError),
            ({**_kwargs_single_valued, "w": -1.0 * u.rad / u.s}, ValueError),
            ({**_kwargs_single_valued, "w": [-1, 2] * u.rad / u.s}, ValueError),
            ({**_kwargs_single_valued, "w": np.ones((2, 2)) * u.rad / u.s}, ValueError),
            ({**_kwargs_single_valued, "w": 5 * u.s}, u.UnitTypeError),
            (
                {**_kwargs_single_valued, "ions": {"not": "a particle"}},
                InvalidParticleError,
            ),
            ({**_kwargs_single_valued, "n_i": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "n_i": 6 * u.m / u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError),
            (
                {**_kwargs_single_valued, "theta": np.ones((2, 2)) * u.rad},
                TypeError,
            ),
            ({**_kwargs_single_valued, "ions": Particle("e-")}, ValueError),
            ({**_kwargs_single_valued, "n_i": [4, 2, 3] * u.m**-3}, ValueError),
            (
                {**_kwargs_single_valued, "n_i": np.ones((2, 2)) * u.m**-3},
                ValueError,
            ),
        ],
    )
    def test_raises(self, kwargs, _error):
        with pytest.raises(_error):
            stix(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({**_kwargs_single_valued, "w": 2 * u.rad / u.s}, {"shape": (4,)}),
            (
                {**_kwargs_single_valued, "w": [10] * u.rad / u.s},
                {"shape": (4,)},
            ),
            (
                {**_kwargs_single_valued, "w": [10, 20, 30] * u.rad / u.s},
                {"shape": (3, 4)},
            ),
            ({**_kwargs_single_valued, "ions": ["He+", "H+"]}, {"shape": (4,)}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": ["He+"],
                    "n_i": [1] * u.m**-3,
                },
                {"shape": (4,)},
            ),
            ({**_kwargs_single_valued, "ions": ["He+", "H+"]}, {"shape": (4,)}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": ["He+", "H+"],
                    "n_i": [1, 2] * u.m**-3,
                },
                {"shape": (4,)},
            ),
            ({**_kwargs_single_valued, "theta": [10, 20, 30]}, {"shape": (3, 4)}),
            (
                {
                    **_kwargs_single_valued,
                    "w": [10, 20] * u.rad / u.s,
                    "theta": 10 * u.rad,
                },
                {"shape": (2, 4)},
            ),
            (
                {
                    **_kwargs_single_valued,
                    "w": [10, 20] * u.rad / u.s,
                    "theta": [0, np.pi / 2, np.pi] * u.rad,
                },
                {"shape": (2, 3, 4)},
            ),
        ],
    )
    def test_return_structure(self, kwargs, expected):
        k = stix(**kwargs)

        assert isinstance(k, u.Quantity)
        assert np.shape(k) == expected["shape"]
        assert k.unit == u.rad / u.m

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            # case taken from Stix figure 1-1
            # Note: ns = [n = k * c / w, ]
            (
                {
                    "theta": 0 * u.rad,
                    "ions": [Particle("p")],
                    "n_i": 1e12 * u.cm**-3,
                    "B": 0.434634 * u.T,
                    "w": 4.16321e4 * u.rad / u.s,
                },
                {
                    "gamma": 1000,
                    "beta": 1000,
                    "mu": 1836,
                    "ns": np.array([31.63146, -31.63146, 31.66306, -31.66306]),
                },
            ),
            #
            # case taken from Stix figure 1-2
            (
                {
                    "theta": 0 * u.rad,
                    "ions": [Particle("p")],
                    "n_i": 1e12 * u.cm**-3,
                    "B": 0.434634 * u.T,
                    "w": (41632 * 10**3) * u.rad / u.s,
                },
                {
                    "gamma": 1000,
                    "beta": 1.1,
                    "mu": 1836,
                    "ns": np.array(
                        [22.39535727, -22.39535727, 6934.82540607, -6934.82540607]
                    ),
                },
            ),
            #
            # case taken from Stix figure 1-3
            (
                {
                    "theta": 0 * u.rad,
                    "ions": [Particle("p")],
                    "n_i": 1e12 * u.cm**-3,
                    "B": 0.434634 * u.T,
                    "w": (124896 * 10**5) * u.rad / u.s,
                },
                {
                    "gamma": 1000,
                    "beta": 1 / 400,
                    "mu": 1836,
                    "ns": np.array(
                        [0.0 + 1.36982834j, 0.0 - 1.36982834j, 2.23009361, -2.23009361]
                    ),
                },
            ),
            #
            # case taken from Stix figure 1-4
            (
                {
                    "theta": 0 * u.rad,
                    "ions": [Particle("p")],
                    "n_i": 1e12 * u.cm**-3,
                    "B": 0.434634 * u.T,
                    "w": (4136 * 10**7) * u.rad / u.s,
                },
                {
                    "gamma": 1000,
                    "beta": 1 / 1300,
                    "mu": 1836,
                    "ns": np.array([0.5880416, -0.5880416, 1.78668573, -1.78668573]),
                },
            ),
            #
            # case taken from Stix figure 1-5
            (
                {
                    "theta": 0 * u.rad,
                    "ions": [Particle("p")],
                    "n_i": 1e12 * u.cm**-3,
                    "B": 0.434634 * u.T,
                    "w": (4136 * 10**7) * u.rad / u.s,
                },
                {
                    "gamma": 1000,
                    "beta": 1 / 1418.5,
                    "mu": 1836,
                    "ns": np.array([0.5880416, -0.5880416, 1.78668573, -1.78668573]),
                },
            ),
        ],
    )
    def test_vals_stix_figs(self, kwargs, expected):
        ion = kwargs["ions"][0]

        mu = ion.mass / Particle("e-").mass
        if not np.isclose(mu, expected["mu"], rtol=9.0e-5):
            pytest.fail(
                "Test setup failure. Check 'kwarg' parameters, given"
                f" values produces a mu of {mu:.2f} but expected "
                f"{expected['mu']:.2f}."
            )

        wpi = plasma_frequency(n=kwargs["n_i"], particle=ion).value
        wci = gyrofrequency(kwargs["B"], particle=ion).value
        gamma = (wpi / wci) ** 2
        if not np.isclose(gamma, expected["gamma"], rtol=0.3e-4):
            pytest.fail(
                "Test setup failure. Check 'kwarg' parameters, given"
                f" values produces a gamma of {gamma:.2f} but expected "
                f"{expected['gamma']:.2f}."
            )

        beta = wci / kwargs["w"].value
        if not np.isclose(beta, expected["beta"], rtol=1):
            pytest.fail(
                "Test setup failure. Check 'kwarg' parameters, given"
                f" values produces a beta of {beta:.3f} but expected "
                f"{expected['beta']:.3f}."
            )

        ks = stix(**kwargs)
        ns = ks.value * c_si_unitless / kwargs["w"].value

        assert np.allclose(ns, expected["ns"])

    @pytest.mark.parametrize(
        "kwargs",
        [
            {
                "ions": ParticleList([Particle("p")]),
                "n_i": [1e12] * u.cm**-3,
                "B": 0.434634 * u.T,
                "w": 4136e7 * u.rad / u.s,
            },
            {
                "ions": ParticleList([Particle("p")]),
                "n_i": [1e12] * u.cm**-3,
                "B": 0.300 * u.T,
                "w": 6e5 * u.rad / u.s,
            },
            {
                "ions": ParticleList([Particle("p")]),
                "n_i": [1e12] * u.cm**-3,
                "B": 0.300 * u.T,
                "w": np.linspace(6e5, 1e9, 10) * u.rad / u.s,
            },
            {
                "ions": ParticleList([Particle("p"), Particle("He+")]),
                "n_i": [
                    0.3 * 1e13,
                    0.7 * 1e13,
                ]
                * u.cm**-3,
                "B": 0.400 * u.T,
                "w": np.linspace(6e5, 1e9, 10) * u.rad / u.s,
            },
        ],
    )
    def test_vals_theta_zero(self, kwargs):
        """
        Test on the known solutions for theta = 0,
        see Stix ch. 1 eqn 37.
        """
        S, P, D = self.spd(**kwargs)
        R = S + D
        L = S - D

        ks = stix(**{**kwargs, "theta": 0 * u.rad})
        ns = ks.value * c_si_unitless / np.tile(kwargs["w"].value, (4, 1)).transpose()

        n_soln1 = np.emath.sqrt(0.5 * (R + L - np.abs(R - L)))  # n^2 = L
        n_soln2 = np.emath.sqrt(0.5 * (R + L + np.abs(R - L)))  # n^2 = R

        assert np.allclose(ns[..., 0], n_soln1)
        assert np.allclose(ns[..., 1], -n_soln1)
        assert np.allclose(ns[..., 2], n_soln2)
        assert np.allclose(ns[..., 3], -n_soln2)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {
                "ions": ParticleList([Particle("p")]),
                "n_i": [1e12] * u.cm**-3,
                "B": 0.434634 * u.T,
                "w": 4136e7 * u.rad / u.s,
            },
            {
                "ions": ParticleList([Particle("p")]),
                "n_i": [1e12] * u.cm**-3,
                "B": 0.300 * u.T,
                "w": 6e5 * u.rad / u.s,
            },
            {
                "ions": ParticleList([Particle("p")]),
                "n_i": [1e12] * u.cm**-3,
                "B": 0.300 * u.T,
                "w": np.linspace(6e5, 1e9, 10) * u.rad / u.s,
            },
            {
                "ions": ParticleList([Particle("p"), Particle("He+")]),
                "n_i": [
                    0.3 * 1e13,
                    0.7 * 1e13,
                ]
                * u.cm**-3,
                "B": 0.400 * u.T,
                "w": np.linspace(6e5, 1e9, 10) * u.rad / u.s,
            },
        ],
    )
    def test_vals_theta_90deg(self, kwargs):
        """
        Test on the known solutions for theta = pi/2,
        see Stix ch. 1 eqn 38.
        """
        S, P, D = self.spd(**kwargs)
        R = S + D
        L = S - D

        ks = stix(**{**kwargs, "theta": 0.5 * np.pi * u.rad})
        ns = ks.value * c_si_unitless / np.tile(kwargs["w"].value, (4, 1)).transpose()

        n_soln1 = np.emath.sqrt(
            (R * L + P * S + np.abs(R * L - P * S)) / (2 * S)
        )  # n^2 = RL / S
        n_soln2 = np.emath.sqrt(
            (R * L + P * S - np.abs(R * L - P * S)) / (2 * S)
        )  # n^2 = P

        assert np.allclose(ns[..., 0], n_soln1)
        assert np.allclose(ns[..., 1], -n_soln1)
        assert np.allclose(ns[..., 2], n_soln2)
        assert np.allclose(ns[..., 3], -n_soln2)
