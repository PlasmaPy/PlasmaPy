from plasmapy.physics.dispersion import two_fluid
import numpy as np

import pytest


tfps_sq_root_table = [
    ((0.6, 1, 0.000545, 0, 1),
     (1.61739, 0.61794, 0.54772)),
    ((1, 2.5, 0.003, 10, 1),
     (3.99630, 1.80755, 1.47895)),
    ((1.5, 2, 0.1, 135, 2),
     (5.84255, 3.19524, 1.06034)),
]


@pytest.mark.parametrize("input_vals, expected", tfps_sq_root_table)
def test_two_fluid_tfps_sq_roots(input_vals, expected):
    sq_roots, _ = two_fluid.tfps(*input_vals)
    assert np.allclose(sq_roots, expected)


tfps_phase_speeds_table = [
    ((0.6, 1, 0.000545, 0, 1),
     (1.61739, 0.61794, 0.54772)),
    ((1, 2.5, 0.003, 10, 1),
     (3.99630, 1.80755, 1.47895)),
    ((1.5, 2, 0.1, 135, 2),
     (2.92127, 1.59762, 0.53017)),
]


@pytest.mark.parametrize("input_vals, expected", tfps_phase_speeds_table)
def test_two_fluid_phase_speeds(input_vals, expected):
    _, phase_speeds = two_fluid.tfps(*input_vals)
    assert np.allclose(phase_speeds, expected)
