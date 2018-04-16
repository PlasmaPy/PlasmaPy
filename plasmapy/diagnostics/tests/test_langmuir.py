# coding=utf-8
"""Tests for Langmuir probe analysis functions."""

import numpy as np
from astropy import units as u
import pytest
import plasmapy.constants as const

from plasmapy.diagnostics.langmuir import (Characteristic,
                                           swept_probe_analysis)

N = 30  # array length of dummy probe characteristic

current_arr = np.random.rand(N) * u.A
bias_arr = np.random.rand(N) * u.V


class Test__characteristic_errors:
    r"""Test the Characteristic class constructor in langmuir.py"""

    bias_infarr = np.append(np.random.rand(N - 1), np.inf) * u.V
    current_infarr = np.append(np.random.rand(N - 1), np.inf) * u.A

    bias_longarr = np.random.rand(N + 1) * u.V
    current_longarr = np.random.rand(N + 1) * u.A

    def test_infinite_I(self):
        r"""Test error upon an infinite current value"""

        with pytest.raises(ValueError):
            Characteristic(bias_arr, self.current_infarr)

    def test_infinite_U(self):
        r"""Test error upon an infinite bias value"""

        with pytest.raises(ValueError):
            Characteristic(self.bias_infarr, current_arr)

    def test_unequal_arrays(self):
        r"""Test errors upon unequal array lengths"""

        with pytest.raises(ValueError):
            Characteristic(self.bias_longarr, current_arr)

        with pytest.raises(ValueError):
            Characteristic(bias_arr, self.current_longarr)


@pytest.fixture
def characteristic():
    r""""Create a dummy characteristic with random values"""
    return Characteristic(bias_arr, current_arr)


@pytest.fixture
def characteristic_simulated():
    r""""Create a simulated probe characteristic (provisional)"""

    T_e_sim = 1 * u.eV
    n_e_sim = 1e18 * u.m**-3
    probe_area_sim = 1 * u.mm**2
    I_e_sim = (n_e_sim * probe_area_sim * const.e *
               np.sqrt(T_e_sim / (2 * np.pi * const.m_e)))

    # Create bias array
    bias_simarr = np.arange(-20, 15, 0.1) * u.V

    # Calculate electron current and limit to electron saturation current
    current_simarr = np.exp(const.e * bias_simarr / T_e_sim) * u.A
    current_simarr[current_simarr > I_e_sim] = I_e_sim

    # Add simulated linear sheath expansion current
    current_simarr[current_simarr == I_e_sim] += (
        bias_simarr[current_simarr == I_e_sim] * 0.0005 * u.A/u.V)

    # Add simulated linear ion collection current
    current_simarr[current_simarr < I_e_sim] += (
        bias_simarr[current_simarr < I_e_sim] * 0.0001 * u.A/u.V)

    return Characteristic(bias_simarr, current_simarr)


@pytest.fixture
def shuffle_characteristic(characteristic):
    r""""Shuffle a given characteristic"""

    _shuffle = sorted(np.arange(len(characteristic.bias)),
                      key=lambda k: np.random.random())
    U_shuffled = characteristic.bias[_shuffle]
    I_shuffled = characteristic.current[_shuffle]
    return Characteristic(U_shuffled, I_shuffled)


class Test__swept_probe_analysis:
    r"""Test the swept_probe_analysis function in langmuir.py"""

    @staticmethod
    def test_nan_area():
        r"""Test error upon NaN area"""

        with pytest.raises(ValueError):
            swept_probe_analysis(characteristic, np.nan * u.m**2, 40 * u.u)

    @staticmethod
    def test_unit_conversion_error():
        r"""Test error upon incorrect probe area unit"""

        with pytest.raises(u.UnitConversionError):
            swept_probe_analysis(characteristic, 1*u.m, 40 * u.u)

    @staticmethod
    def test_negative_area():
        r"""Test error upon negative probe area"""

        with pytest.raises(ValueError):
            swept_probe_analysis(characteristic, -1 * u.m**2, 40 * u.u)

    @staticmethod
    def test_ion_mass_unit():
        r"""Test equality of float and a.m.u. ion mass"""

        sim = characteristic_simulated()

        sim_result1 = swept_probe_analysis(
            sim, 4*u.m**2, 40)

        sim_result2 = swept_probe_analysis(
            sim, 4*u.m**2, 40 * u.u)

        errStr = (f"`swept_probe_analysis` should accept both floats and "
                  f"a.m.u. Quantities as atomic gas mass input.")
        assert sim_result1 == sim_result2, errStr

    @staticmethod
    def test_ordering_invariance():
        r"""Test invariance to ordering of the bias and current values"""

        sim = characteristic_simulated()

        sim_result = swept_probe_analysis(
            sim, 4*u.m**2, 40 * u.u)

        sim_result_shuffled = swept_probe_analysis(
            shuffle_characteristic(sim), 4*u.m**2, 40 * u.u)

        errStr = (f"Analysis should be invariant to the ordering of the "
                  f"input data.")
        assert sim_result == sim_result_shuffled, errStr
