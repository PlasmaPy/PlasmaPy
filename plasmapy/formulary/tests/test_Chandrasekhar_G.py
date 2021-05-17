import hypothesis
import numpy as np

from hypothesis import example, given
from hypothesis import strategies as st

from plasmapy.formulary.mathematics import Chandrasekhar_G


@given(
    x=st.floats(
        allow_nan=False,
        allow_infinity=False,
    )
)
@example(x=np.finfo(np.float64).eps)
@example(x=np.finfo(np.float64).max)
@example(x=np.finfo(np.float64).min)
@example(x=3.4694469519536216e-17)
@example(x=0)
@example(x=3.761264e-20)
def test_Chandrasekhar_with_hypothesis(x):
    with np.errstate():
        result = Chandrasekhar_G(x)
    assert abs(result) < 0.21399915915288345  # maximum bound found via scipy optimize
    assert np.isfinite(result)

    symmetric_result = Chandrasekhar_G(-x)  # antisymmetric function
    assert result == -symmetric_result


def test_Chandrasekhar_regression(num_regression):
    x = np.linspace(-5, 5, 100)
    y = Chandrasekhar_G(x)
    num_regression.check({"x": x, "y": y})


def test_limiting_behavior_small_x():
    x = np.logspace(-7, -4, 100)
    y = Chandrasekhar_G(x)
    limiting_behavior = 2 * x / (3 * np.sqrt(np.pi))
    # useful debug check
    # import matplotlib.pyplot as plt
    # plt.semilogx(x, y - limiting_behavior)
    # plt.show()
    np.testing.assert_allclose(y, limiting_behavior, rtol=0, atol=1e-8)


def test_limiting_behavior_high_x():
    x = np.logspace(6, 90, 100)
    y = Chandrasekhar_G(x)
    limiting_behavior = x ** -2 / 2
    np.testing.assert_allclose(y, limiting_behavior)
