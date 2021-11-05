from astropy import units as u
import numpy as np
import pytest
from stix_ import stix

class StixTest:
	_kwargs_single_valued = {
		"B": 8.3e-9 * u.T,
		"k": 0.001* u.rad / u.m,
		"ions": ['e-','H+'],
		"omega_ions": [4.0e5, 2.0e5] * u.rad / u.s,
		"theta": 30 * u.deg,
	}

	@pytest.mark.parametrize(
		"kwargs, _error",
		[
			({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
			({**_kwargs_single_valued, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
			({**_kwargs_single_valued, "B": -1 * u.T}, ValueError),
            		({**_kwargs_single_valued, "B": 5 * u.m}, u.UnitTypeError),
            		({**_kwargs_single_valued, "k": 0 * u.rad / u.m}, ValueError),
            		({**_kwargs_single_valued, "k": -1.0 * u.rad / u.m}, ValueError),
            		({**_kwargs_single_valued, "k": 5 * u.s}, u.UnitTypeError),
            		({**_kwargs_single_valued, "ions": {"not": "a particle"}}, TypeError),
           		({**_kwargs_single_valued, "ions": "e-"}, ValueError),
			({**_kwargs_single_valued, "omega_ions": "wrong type"}, TypeError),
			({**_kwargs_single_valued, "omega_ions": 6 * u.eV}, u.UnitTypeError),
			({**_kwargs_single_valued, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            		({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError)
		]
	)
	def test_raises(self, kwargs, _error):
		with pytest.raises(_error):
            		stix(**kwargs)

       
print('Out')
