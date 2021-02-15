import numpy as np

from plasmapy.formulary.mathematics import Chandrasekhar_G


def test_Chandrasekhar(num_regression):
    x = np.linspace(-5, 5, 100)
    y = Chandrasekhar_G(x)
    num_regression.check({"x": x, "y": y})
