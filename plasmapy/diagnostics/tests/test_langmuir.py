# coding=utf-8
"""Tests for Langmuir probe analysis functions."""

import numpy as np
from astropy import units as u
import pytest
import plasmapy.constants as const

from plasmapy.diagnostics.langmuir import (Characteristic,
                                           swept_probe_analysis)

N = 30 # length of test characteristic

I_arr = np.random.rand(N) * u.A
I_infarr = np.append(np.random.rand(N - 1), np.inf) * u.A

U_arr = np.random.rand(N) * u.V
U_infarr = np.append(np.random.rand(N - 1), np.inf) * u.V

def test_characteristic():
    r"""Test the Characteristic class constructor in langmuir.py"""
        
    with pytest.raises(ValueError):
        Characteristic(U_arr, I_infarr)
        
    with pytest.raises(ValueError):
        Characteristic(U_infarr, I_arr)

def test_swept_probe_analysis():
    r"""Test the swept_probe_analysis function in langmuir.py"""
    
    with pytest.raises(ValueError):
        swept_probe_analysis(Characteristic(U_arr, I_arr), np.nan * u.m**2, 40)
    
    with pytest.raises(u.UnitConversionError):
        swept_probe_analysis(Characteristic(U_arr, I_arr), 1*u.m, 40)
    
    with pytest.raises(ValueError):
        swept_probe_analysis(Characteristic(U_arr, I_arr), -1 * u.m**2, 40)
        
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
    
    _shuffle = sorted(np.arange(len(U_simarr)), key=lambda k: np.random.random())
    U_simarr_shuffled = U_simarr[_shuffle]
    I_simarr_shuffled = I_simarr[_shuffle]
    
    sim_result = swept_probe_analysis(Characteristic(U_simarr, I_simarr), 
                                      4*u.m**2, 40)
    
    sim_result_shuffled = swept_probe_analysis(Characteristic(U_simarr_shuffled, 
                                               I_simarr_shuffled), 
                                               4*u.m**2, 40)
    
    errStr = (f"Analysis should be invariant to the ordering of the input "
              f"data")
    assert sim_result == sim_result_shuffled, errStr