"""Test that objects in `plasmapy.particles` can be pickled."""

import astropy.units as u
import os
import pickle
import pytest

from plasmapy.particles.ionization_state import IonicLevel, IonizationState
from plasmapy.particles.ionization_state_collection import IonizationStateCollection
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
)


class TestPickling:
    """
    Test that different objects in `plasmapy.particles` can be pickled.
    """

    xfail = pytest.mark.xfail(reason="see issue #1011")

    @pytest.mark.parametrize(
        "instance",
        [
            CustomParticle(mass=1 * u.kg, charge=1 * u.C),
            DimensionlessParticle(mass=5, charge=5),
            pytest.param(Particle("p+"), marks=xfail),
            pytest.param(IonicLevel("p+", 0.1, 1e9 * u.m ** -3), marks=xfail),
            pytest.param(IonizationState("H", [0.5, 0.5]), marks=xfail),
            pytest.param(IonizationStateCollection({"H": [0.5, 0.5]}), marks=xfail),
        ],
    )
    def test_pickling(self, instance, tmp_path):
        """
        Test that different objects contained within `plasmapy.particles`
        can be pickled and unpickled.
        """
        filename = tmp_path / "pickled_particles.p"
        with open(filename, "wb") as pickle_file:
            pickle.dump(instance, pickle_file)

        with open(filename, "rb") as pickle_file:
            loaded_particle = pickle.load(pickle_file)

        assert str(instance) == str(loaded_particle)
