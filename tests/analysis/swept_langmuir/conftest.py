import numpy as np
import pytest


@pytest.fixture(scope="module")
def simple_voltage():
    return np.linspace(-20.0, 20.0, 200)


@pytest.fixture(scope="module")
def simple_current(simple_voltage):
    return simple_voltage.copy()


def arch_current(voltage, *, radius, centroid, positive_arch=True):
    # radius: float
    #     radius of the arch/circle to be generated
    #
    # centroid: 2 element array
    #     center (x, y) location of the circle
    #
    # positive_arch: bool
    #     if True draw the positive (upper) arch
    #     if False draw the negative (lower) arch
    #

    sign = 1.0 if positive_arch else -1.0
    return sign * np.sqrt(
        radius ** 2 - (voltage - centroid[0]) ** 2
    ) + centroid[1]
