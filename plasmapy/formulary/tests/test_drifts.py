import astropy.units as u
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary import drifts


class Test_diamagnetic_drift:
    q = 1 * u.C
    n = 1 / u.m**3

    def test_isothermal_plasma(self):
        r"""
        Cylindrical & isothermal plasma in Fig 3.4 of Chen.
        """
        T = 10 * u.eV
        B = u.Quantity([0, 0, 100], unit=u.nT)
        # Lambda = n' / n = 1 / m
        dp = T * u.Quantity([0, -self.n.value, 0], unit=self.n.unit / u.m)

        result = drifts.diamagnetic_drift(dp, B, self.n, self.q)
        assert_quantity_allclose(
            np.linalg.norm(result), T / self.q / np.linalg.norm(B) / u.m
        )

    def test_diamagnetic_drift_1d_arrays(self):
        dp = u.Quantity([1, 0, 0], unit=u.N / u.m**3)
        B = u.Quantity([0, -1, 0], unit=u.T)
        result = drifts.diamagnetic_drift(2 * dp, 3 * B, self.n, self.q)
        assert_quantity_allclose(result, (2 / 3) * u.Quantity([0, 0, 1], u.m / u.s))

    def test_diamagnetic_drift_2d_array(self):
        dp = u.Quantity([[1, 0, 0], [1, 0, 0], [1, 0, 0]], unit=u.N / u.m**3)
        B = -u.Quantity([[0, 1, 0], [0, 1, 0], [0, 1, 0]], unit=u.T)

        result = drifts.diamagnetic_drift(2 * dp, 3 * B, self.n, self.q)
        assert_quantity_allclose(
            result,
            (2 / 3) * u.Quantity([[0, 0, 1], [0, 0, 1], [0, 0, 1]], unit=u.m / u.s),
        )

    def test_diamagnetic_drift_3d_array(self):
        dp = u.Quantity([[[1, 0, 0]]], unit=u.N / u.m**3)
        B = u.Quantity([[[0, -1, 0]]], unit=u.T)

        result = drifts.diamagnetic_drift(2 * dp, 3 * B, self.n, self.q)
        assert_quantity_allclose(
            result, (2 / 3) * u.Quantity([[[0, 0, 1]]], unit=u.m / u.s)
        )

    def test_nonsensical_units(self):
        dp = u.Quantity([[1, 0, 0], [1, 0, 0], [1, 0, 0]], unit=u.C)
        B = u.Quantity([[0, 1, 0], [0, 1, 0], [0, 1, 0]], unit=u.kg)
        q = 1 * u.C
        n = 1 / u.m**3

        with pytest.raises(u.UnitTypeError):
            drifts.diamagnetic_drift(dp, B, n, q)

    def test_alias(self):
        assert drifts.vd_ is drifts.diamagnetic_drift


class Test_ExB_drift:
    def test_E_x_B_1d_arrays(self):
        E = u.Quantity([1, 0, 0], unit=u.V / u.m)
        B = u.Quantity([0, 1, 0], unit=u.T)
        result = drifts.ExB_drift(2 * E, 3 * B)
        assert_quantity_allclose(result, (2 / 3) * u.Quantity([0, 0, 1], u.m / u.s))

    def test_ExB_2d_array(self):
        E = u.Quantity([[1, 0, 0], [1, 0, 0], [1, 0, 0]], unit=u.V / u.m)
        B = u.Quantity([[0, 1, 0], [0, 1, 0], [0, 1, 0]], unit=u.T)

        result = drifts.ExB_drift(2 * E, 3 * B)
        assert_quantity_allclose(
            result,
            (2 / 3) * u.Quantity([[0, 0, 1], [0, 0, 1], [0, 0, 1]], unit=u.m / u.s),
        )

    def test_nonsensical_units(self):
        E = u.Quantity([[1, 0, 0], [1, 0, 0], [1, 0, 0]], unit=u.mm)
        B = u.Quantity([[0, 1, 0], [0, 1, 0], [0, 1, 0]], unit=u.kg)

        with pytest.raises(u.UnitTypeError):
            drifts.ExB_drift(E, B)

    def test_ExB_3d_array(self):
        E = u.Quantity([[[1, 0, 0]]], unit=u.V / u.m)
        B = u.Quantity([[[0, 1, 0]]], unit=u.T)

        result = drifts.ExB_drift(2 * E, 3 * B)
        assert_quantity_allclose(
            result, (2 / 3) * u.Quantity([[[0, 0, 1]]], unit=u.m / u.s)
        )

    def test_alias(self):
        assert drifts.veb_ is drifts.ExB_drift


class Test_force_drift:
    def test_force_x_B_1d_arrays(self):
        F = u.Quantity([1, 0, 0], unit=u.N)
        B = u.Quantity([0, 1, 0], unit=u.T)
        q = 1 * u.C
        result = drifts.force_drift(2 * F, 3 * B, q)
        assert_quantity_allclose(result, (2 / 3) * u.Quantity([0, 0, 1], u.m / u.s))

    def test_force_B_2d_array(self):
        F = u.Quantity([[1, 0, 0], [1, 0, 0], [1, 0, 0]], unit=u.N)
        B = u.Quantity([[0, 1, 0], [0, 1, 0], [0, 1, 0]], unit=u.T)
        q = 1 * u.C

        result = drifts.force_drift(2 * F, 3 * B, q)
        assert_quantity_allclose(
            result,
            (2 / 3) * u.Quantity([[0, 0, 1], [0, 0, 1], [0, 0, 1]], unit=u.m / u.s),
        )

    def test_force_B_3d_array(self):
        F = u.Quantity([[[1, 0, 0]]], unit=u.N)
        B = u.Quantity([[[0, 1, 0]]], unit=u.T)
        q = 1 * u.C

        result = drifts.force_drift(2 * F, 3 * B, q)
        assert_quantity_allclose(
            result, (2 / 3) * u.Quantity([[[0, 0, 1]]], unit=u.m / u.s)
        )

    def test_nonsensical_units(self):
        F = u.Quantity([[1, 0, 0], [1, 0, 0], [1, 0, 0]], unit=u.mm)
        B = u.Quantity([[0, 1, 0], [0, 1, 0], [0, 1, 0]], unit=u.kg)
        q = 1 * u.C

        with pytest.raises(u.UnitTypeError):
            drifts.force_drift(F, B, q)

    def test_alias(self):
        assert drifts.vfd_ is drifts.force_drift
