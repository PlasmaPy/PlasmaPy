# coding=utf-8
"""Tests for Langmuir probe analysis functions."""

import numpy as np
from astropy import units as u
import pytest
import plasmapy.constants as const

from plasmapy.diagnostics import langmuir

np.random.seed(42)
N = 30  # array length of dummy probe characteristic

current_arr = np.random.rand(N) * u.A
bias_arr = np.random.rand(N) * u.V


class Test__fitting_functions:
    x = 20
    x0 = 10
    y0 = 5
    c0 = 2
    T0 = 3
    Delta_T = 15

    def test_fit_func_lin(self):
        value = langmuir._fit_func_lin(self.x, self.x0, self.y0, self.c0)
        assert value == 25

    def test_fit_func_lin_inverse(self):
        value = langmuir._fit_func_lin_inverse(self.x, self.x0, self.y0, self.T0)
        assert np.allclose(value, 8.33333)

    def test_fit_func_double_lin_inverse(self):
        value = langmuir._fit_func_double_lin_inverse(self.x, self.x0, self.y0,
                                                      self.T0, self.Delta_T)
        assert value == 8


class Test__characteristic_errors:
    r"""Test the Characteristic class constructor in langmuir.py"""

    bias_2darr = np.array((np.random.rand(N),
                           np.random.rand(N))) * u.V
    current_2darr = np.array((np.random.rand(N),
                              np.random.rand(N))) * u.A

    bias_infarr = np.append(np.random.rand(N - 1), np.inf) * u.V
    current_infarr = np.append(np.random.rand(N - 1), np.inf) * u.A

    bias_longarr = np.random.rand(N + 1) * u.V
    current_longarr = np.random.rand(N + 1) * u.A

    current_arr2 = np.random.rand(N) * u.A

    def test_invalid_dimensions(self):
        # Checks if `Characteristic.get_unique_bias` runs during
        # initialization.
        langmuir.Characteristic(self.bias_2darr, self.current_2darr)

    def test_infinite_I(self):
        r"""Test error upon an infinite current value"""

        with pytest.raises(ValueError):
            langmuir.Characteristic(bias_arr, self.current_infarr)

    def test_infinite_U(self):
        r"""Test error upon an infinite bias value"""

        with pytest.raises(ValueError):
            langmuir.Characteristic(self.bias_infarr, current_arr)

    def test_unequal_arrays(self):
        r"""Test errors upon unequal array lengths"""

        with pytest.raises(ValueError):
            langmuir.Characteristic(self.bias_longarr, current_arr)

        with pytest.raises(ValueError):
            langmuir.Characteristic(bias_arr, self.current_longarr)

    def test_addition(self):
        r"""Test addition of characteristic objects"""

        a = langmuir.Characteristic(bias_arr, current_arr)

        b = langmuir.Characteristic(bias_arr, self.current_arr2)

        ab_sum = a + b

        errStr = (f"Addition of characteristic objects is not behaving as it "
                  f"should.")
        assert (a.current + b.current == ab_sum.current).all(), errStr

    def test_subtraction(self):
        r"""Test addition of characteristic objects"""

        a = langmuir.Characteristic(bias_arr, current_arr)

        b = langmuir.Characteristic(bias_arr, self.current_arr2)

        ab_sub = a - b

        errStr = (f"Subtraction of characteristic objects is not behaving as "
                  f"it should.")
        assert (a.current - b.current == ab_sub.current).all(), errStr


@pytest.fixture
def characteristic():
    r""""Create a dummy characteristic with random values"""
    return langmuir.Characteristic(bias_arr, current_arr)


@pytest.fixture
def characteristic_simulated():
    r""""Create a simulated probe characteristic (provisional)"""

    T_e_sim = 1 * u.eV
    n_e_sim = 1e18 * u.m**-3
    probe_area_sim = 1 * u.cm**2
    I_es_sim = (n_e_sim * probe_area_sim * const.e *
                np.sqrt(T_e_sim / (2 * np.pi * const.m_e)))

    # Create bias array
    bias_simarr = np.arange(-20, 15, 0.1) * u.V

    # Calculate electron current and limit to electron saturation current
    current_simarr = np.exp(const.e * bias_simarr / T_e_sim) * u.A
    current_simarr[current_simarr > I_es_sim] = I_es_sim

    # Add simulated linear sheath expansion current
    current_simarr[current_simarr == I_es_sim] += (
        bias_simarr[current_simarr == I_es_sim] * 5e-4 * u.A/u.V)

    # Add simulated linear ion collection current
    current_simarr[current_simarr < I_es_sim] += (
        bias_simarr[current_simarr < I_es_sim] * 1e-4 * u.A/u.V)

    return langmuir.Characteristic(bias_simarr, current_simarr)


@pytest.fixture
def shuffle_characteristic(characteristic):
    r""""Shuffle a given characteristic"""

    _shuffle = sorted(np.arange(len(characteristic.bias)),
                      key=lambda k: np.random.random())
    U_shuffled = characteristic.bias[_shuffle]
    I_shuffled = characteristic.current[_shuffle]
    return langmuir.Characteristic(U_shuffled, I_shuffled)


class DryCharacteristic(langmuir.Characteristic):
    def __init__(self, bias, current):
        self.bias = bias
        self.current = current


class Test__Characteristic_methods:
    r"""."""
    bias_2darr = np.array((np.random.rand(N),
                           np.random.rand(N))) * u.V
    current_2darr = np.array((np.random.rand(N),
                              np.random.rand(N))) * u.A

    bias_4length_arr = np.array(np.random.rand(N-1)) * u.V
    current_5length_arr = np.array(np.random.rand(N)) * u.A

    bias_duplicates_arr = np.array((1,2) * int(N/2))

    def test_invalid_bias_dimensions(self):
        with pytest.raises(ValueError):
            char = DryCharacteristic(self.bias_2darr,
                                     current_arr)
            char.check_validity()

    def test_invalid_current_dimensions(self):
        with pytest.raises(ValueError):
            char = DryCharacteristic(bias_arr,
                                     self.current_2darr)
            char.check_validity()

    def test_bias_and_current_length_mismatch(self):
        with pytest.raises(ValueError):
            char = DryCharacteristic(self.bias_4length_arr,
                                     self.current_5length_arr)
            char.check_validity()

    def test_duplicate_bias_values(self):
        with pytest.raises(ValueError):
            char= DryCharacteristic(self.bias_duplicates_arr,
                                    current_arr)
            char.check_validity()

    @staticmethod
    def test_inplace_unique_bias():
        char = DryCharacteristic(bias_arr, current_arr)
        new_char = char.get_unique_bias(inplace=False)
        assert char != new_char
        assert isinstance(new_char, langmuir.Characteristic)

    @staticmethod
    def test_getpadded_limit():
        char = characteristic()
        limits = char.get_padded_limit(0.1)
        log_limits = char.get_padded_limit(0.1, log=True)
        assert np.allclose(limits.to(u.A).value,
                           np.array((-0.07434804, 1.06484239)))
        assert np.allclose(log_limits.to(u.A).value,
                           np.array((0.014003, 1.42577333)))


class Test__swept_probe_analysis:
    r"""Test the swept_probe_analysis function in langmuir.py"""

    @staticmethod
    def test_nan_area():
        r"""Test error upon NaN area"""

        with pytest.raises(ValueError):
            langmuir.swept_probe_analysis(characteristic, np.nan * u.cm**2, 40 * u.u)

    @staticmethod
    def test_unit_conversion_error():
        r"""Test error upon incorrect probe area unit"""

        with pytest.raises(u.UnitConversionError):
            langmuir.swept_probe_analysis(characteristic, 1 * u.cm, 40 * u.u)

    @staticmethod
    def test_negative_area():
        r"""Test error upon negative probe area"""

        with pytest.raises(ValueError):
            langmuir.swept_probe_analysis(characteristic, -1 * u.cm**2, 40 * u.u)

    @staticmethod
    def test_ion_mass_unit():
        r"""Test equality of float and a.m.u. ion mass"""

        sim = characteristic_simulated()

        sim_result1 = langmuir.swept_probe_analysis(
            sim, 1 * u.cm**2, 40)

        sim_result2 = langmuir.swept_probe_analysis(
            sim, 1 * u.cm**2, 40 * u.u)

        errStr = (f"`swept_probe_analysis` should accept both floats and "
                  f"a.m.u. Quantities as atomic gas mass input.")
        for key in sim_result1:
            assert (sim_result1[key] == sim_result2[key]).all(), errStr

    @staticmethod
    @pytest.mark.parametrize("bimaxwellian", [True, False])
    def test_ordering_invariance(bimaxwellian):
        r"""Test invariance to ordering of the bias and current values"""

        sim = characteristic_simulated()

        sim_result = langmuir.swept_probe_analysis(
            sim,
            1 * u.cm**2,
            40 * u.u,
            bimaxwellian=bimaxwellian)

        sim_result_shuffled = langmuir.swept_probe_analysis(
            shuffle_characteristic(sim),
            1 * u.cm**2,
            40 * u.u,
            bimaxwellian=bimaxwellian)

        errStr = (f"Analysis should be invariant to the ordering of the "
                  f"input data.")
        for key in sim_result:
            assert (sim_result[key] == sim_result_shuffled[key]).all(), errStr
