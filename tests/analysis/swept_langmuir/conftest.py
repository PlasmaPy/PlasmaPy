import numpy as np
import pytest


@pytest.fixture(scope="module")
def simple_voltage():
    return np.linspace(-20.0, 20.0, 100)


@pytest.fixture(scope="module")
def simple_current():
    return np.linspace(-20.0, 20.0, 100)
