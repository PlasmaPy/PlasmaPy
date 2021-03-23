import hypothesis
import numpy as np

from hypothesis import given, settings
from hypothesis import strategies as st

from plasmapy.formulary.mathematics import Chandrasekhar_G


def test_Chandrasekhar(num_regression):
    x = np.linspace(-5, 5, 100)
    y = Chandrasekhar_G(x)
    num_regression.check({"x": x, "y": y})


@given(
    x=st.floats(
        min_value=6e-150,
        max_value=1e90,
        exclude_min=True,
        allow_nan=False,
        allow_infinity=False,
    )
)
@settings(max_examples=500)
def test_Chandrasekhar_with_hypothesis(x):
    result = Chandrasekhar_G(x)
    assert x > 0
    assert np.isfinite(x)
