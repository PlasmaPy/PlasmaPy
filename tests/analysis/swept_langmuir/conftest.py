"""
pytest fixtures for all swept_langmuir tests.
"""

from pathlib import Path

import numpy as np
import pytest

from plasmapy.analysis.swept_langmuir import sort_sweep_arrays

PACE_DATA_FILE = (Path(__file__).parent / "Pace2015.npy").resolve()
PACE_DATA = sort_sweep_arrays(*np.load(PACE_DATA_FILE))


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
    return sign * np.sqrt(radius**2 - (voltage - centroid[0]) ** 2) + centroid[1]


@pytest.fixture(scope="module")
def simple_current_two_arches(simple_voltage):
    voltage = simple_voltage
    current = np.full_like(voltage, np.nan)

    # inflection_point = location where the two arches transition from
    #     one to the other
    # slope_angle = angle of the slope line at the two arch inflection
    #     point
    # bisect_line_angle = the angle of the line that bisects the
    #     centroids of the two arch circles
    #
    inflection_point = np.array([0.0, 0.0])
    slope_angle = np.deg2rad(65.0)
    bisect_line_angle = slope_angle - np.deg2rad(90)

    centroid_1 = np.array([1.1 * np.min(voltage), np.nan])
    centroid_2 = np.array([1.1 * np.max(voltage), np.nan])

    radius_1 = np.abs(centroid_1[0]) / np.cos(bisect_line_angle)
    radius_2 = np.abs(centroid_2[0]) / np.cos(bisect_line_angle)

    centroid_1[1] = -radius_1 * np.sin(bisect_line_angle) + inflection_point[1]
    centroid_2[1] = radius_2 * np.sin(bisect_line_angle) + inflection_point[1]

    # generate first (left) arch
    mask = voltage <= 0
    current[mask] = arch_current(
        voltage[mask],
        centroid=centroid_1,
        radius=radius_1,
        positive_arch=False,
    )

    # generate second (right) arch
    mask = voltage > 0
    current[mask] = arch_current(
        voltage[mask],
        centroid=centroid_2,
        radius=radius_2,
        positive_arch=True,
    )

    return current


@pytest.fixture(scope="module")
def pace_voltage():
    return PACE_DATA[0]


@pytest.fixture(scope="module")
def pace_current():
    return PACE_DATA[1]
