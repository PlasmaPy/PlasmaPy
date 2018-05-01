import pytest
import numpy as np
from ..ionization_states import IonizationState


def test_instance():

    elem = 'Li'
    dist = np.array([0.8, 0.15, 0.05, 0.0])

    ionization_state = IonizationState(
        element=elem,
        distribution=dist,
    )

    # TODO: Expand tests and make them non-serial

    assert ionization_state[0] == dist[0]

    assert np.allclose(ionization_state[0:3], dist[0:3])

    assert np.allclose(ionization_state[0:3:2], dist[0:3:2])