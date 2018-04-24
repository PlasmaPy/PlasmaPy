import pytest
from ..ionization_states import IonizationState


def test_instance():
    ionization_state = IonizationState(
        element='Li',
        distribution=[0.8, 0.15, 0.05, 0.0],
    )

    assert ionization_state[0] == 0.8
