# coding=utf-8
"""Tests for Langmuir probe analysis functions."""

import numpy as np
from astropy import units as u
import pytest
import plasmapy.constants as const

from plasmapy.diagnostics.langmuir import (Characteristic,
                                           swept_probe_analysis)

N = 30  # length of test characteristic

I_arr = np.random.rand(N) * u.A
U_arr = np.random.rand(N) * u.V
class Test__characteristic_errors:
    r"""Test the Characteristic class constructor in langmuir.py"""
    N = 30  # length of test characteristic
    U_infarr = np.append(np.random.rand(N - 1), np.inf) * u.V
    I_infarr = np.append(np.random.rand(N - 1), np.inf) * u.A

    def test_infinite_I(self):
        with pytest.raises(ValueError):
            Characteristic(U_arr, self.I_infarr)

    def test_infinite_U(self):
        with pytest.raises(ValueError):
            Characteristic(self.U_infarr, I_arr)


@pytest.fixture
def characteristic():
    return Characteristic(U_arr, I_arr)

@pytest.fixture
def characteristic_simulated():
    r""""Simulated characteristic check below (unfinished)"""
    T_e_sim = 1 * u.eV
    n_e_sim = 10**18 * u.m**-3
    probe_area_sim = 1 * u.mm**2
    I_e_sim = n_e_sim * probe_area_sim * const.e * \
        np.sqrt(T_e_sim / (2 * np.pi * const.m_e))

    U_simarr = np.arange(-20, 15, 0.1) * u.V
    I_simarr = np.exp(const.e * U_simarr / T_e_sim) * u.A
    I_simarr[I_simarr > I_e_sim] = I_e_sim
    I_simarr[I_simarr < I_e_sim] += U_simarr[I_simarr < I_e_sim] * \
        0.0001 * u.A/u.V
    I_simarr[I_simarr == I_e_sim] += U_simarr[I_simarr == I_e_sim] * \
        0.0005 * u.A/u.V
    return Characteristic(U_simarr, I_simarr)

@pytest.fixture
def characteristic_simulated_shuffle(characteristic_simulated):
    _shuffle = sorted(np.arange(len(characteristic_simulated.bias)),
                      key=lambda k: np.random.random())
    U_simarr_shuffled = characteristic_simulated.bias[_shuffle]
    I_simarr_shuffled = characteristic_simulated.current[_shuffle]
    return Characteristic(U_simarr_shuffled, I_simarr_shuffled)

class Test__swept_probe_analysis:
    r"""Test the swept_probe_analysis function in langmuir.py"""

    @staticmethod
    def test_nan_area(characteristic):
        with pytest.raises(ValueError):
            swept_probe_analysis(characteristic, np.nan * u.m**2, 40)

    @staticmethod
    def test_unit_conversion_error(characteristic):
        with pytest.raises(u.UnitConversionError):
            swept_probe_analysis(characteristic, 1*u.m, 40)

    @staticmethod
    def test_negative_surface(characteristic):
        with pytest.raises(ValueError):
            swept_probe_analysis(characteristic, -1 * u.m**2, 40)

    @staticmethod
    def test_ordering_invariance(characteristic_simulated, characteristic_simulated_shuffle):
        sim_result = swept_probe_analysis(characteristic_simulated, 4*u.m**2, 40)

        sim_result_shuffled = swept_probe_analysis(characteristic_simulated_shuffle, 4*u.m**2, 40)

        errStr = "Analysis should be invariant to the ordering of the input data."
        assert sim_result == sim_result_shuffled, errStr
