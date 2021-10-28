from astropy import units as u

from stix_ import stix

import pytest

B = 8.3e-9 * u.T
k_int = 0.001 * u.rad / u.m
k_arr = [0.001, 0.002] * u.rad / u.m
ions_int = "e-"
ions_arr = ["e-", "H+"]
omega_ions_int = 4.0e5 * u.rad / u.s
omega_ions_arr = [4.0e5, 2.0e5] * u.rad / u.s
theta = 30 * u.deg

inputs = {}
eops = {}

inputs = {
    "B": B,
    "k": k_int,
    "ions": ions_int,
    "omega_ions": omega_int,
    "theta": theta,
}

eops = {

}

def test_confrim(inputs, eops):
    assert stix(**inputs) == eops

inputs = {
    "B": B,
    "k": k_int,
    "ions": ions_arr,
    "omega_ions": omega_arr,
    "theta": theta,
}

eops = {

}

def test_confirm(inputs, eops):
    assert stix(**inputs) == eops

inputs = {
    "B": B,
    "k": k_arr,
    "ions": ions_int,
    "omega_ions": omega_int,
    "theta": theta,
}

eops = {

}

def test_confirm(inputs, eops):
    assert stix(**inputs) == eops

inputs = {
    "B": B,
    "k": k_arr,
    "ions": ions_arr,
    "omega_ions": omega_arr,
    "theta": theta,
}

eops = {

}

def test_confirm(inputs, eops):
    assert stix(**inputs) == eops

inputs = {
    "B": B,
    "k": k_arr,
    "ions": ions_arr,
    "omega_ions": omega_int,
    "theta": theta,
}
    
def test_flag(inputs):
    with pytest.raises(SystemExit):
        stix(**inputs)
