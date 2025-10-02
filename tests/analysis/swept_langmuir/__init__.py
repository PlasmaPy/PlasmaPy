import numpy as np
import pytest


class SLAFixtures:

    @pytest.fixture(scope="class")
    def simple_voltage(self):
        return np.linspace(-20.0, 20.0, 100)

    @pytest.fixture(scope="class")
    def simple_current(self):
        return np.linspace(-20.0, 20.0, 100)
