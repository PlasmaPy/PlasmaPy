"""Tests for the UncertaintyQuantity class"""

import numpy as np
import astropy.units as u

from plasmapy.uncertainty import UncertaintyQuantity

def test_UncertaintyQuantity():

    base0, unc0 = 30, 5
    base1, unc1 = 10, 1.5

    value0 = UncertaintyQuantity(base0 * u.m, unc0 * u.m)
    value1 = UncertaintyQuantity(base1 * u.m, unc1 * u.m)

    assert value0 == UncertaintyQuantity(base0 * u.m, unc0 * u.m)

    assert value0 != value1

    assert value0 / u.m == UncertaintyQuantity(base0, unc0)

    assert 2 * value0 == UncertaintyQuantity(2 * base0 * u.m, 2 * unc0 * u.m)

    assert value0 > value1

    assert value0 == UncertaintyQuantity(base0, unc0) * u.m

    assert np.sqrt(value0) == value0 ** 0.5

    assert np.power(value0, 2) == value0 ** 2

    assert np.power(2, value0 / u.m) == 2 ** (value0 / u.m)

    assert np.add(value0, value1) == value0 + value1

    assert np.add(value0, value1) == np.add(value1, value0)

    assert np.subtract(value0, value1) == value0 - value1

    assert np.multiply(value0, value1) == value0 * value1

    assert np.multiply(value0, value1) == np.multiply(value1, value0)

    assert np.multiply(value0, 1) == value0

    assert np.multiply(1, value1) == value1

    assert np.true_divide(value0, value1) == value0 / value1

    assert np.square(value0) == value0 ** 2

    value0_converted = value0.to(u.imperial.inch)

    assert value0_converted.unit == u.imperial.inch

    assert value0_converted.uncertainty.unit == u.imperial.inch