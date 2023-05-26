"""Tests the Fermi_integral function in the mathematics module."""

# This file contains experimental usage of unicode characters.

import numpy as np
import pytest

from plasmapy.formulary.mathematics import Fermi_integral


class Test_Fermi_integral:
    @classmethod
    def setup_class(cls):
        """Initialize parameters for tests."""
        cls.arg1 = 3.889780
        cls.True1 = 6.272518847136373 - 8.673617379884035e-19j
        cls.argFail1 = 3.889781
        cls.order1 = 0.5
        cls.args = np.array([0.5, 1, 2])
        cls.Trues = np.array(
            [
                (1.1173314873128224 - 0j),
                (1.5756407761513003 - 0j),
                (2.8237212774015843 - 2.6020852139652106e-18j),
            ]
        )

    def test_known1(self):
        """
        Test Fermi_integral for expected value.
        """
        methodVal = Fermi_integral(self.arg1, self.order1)
        testTrue = np.isclose(methodVal, self.True1, rtol=1e-16, atol=0.0)
        errStr = f"Fermi integral value should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Test if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 + 1e-15
        methodVal = Fermi_integral(self.arg1, self.order1)
        testTrue = not np.isclose(methodVal, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Fermi integral value test gives {methodVal} and should "
            f"not be equal to {fail1}."
        )
        assert testTrue, errStr

    def test_array(self):
        """Test Fermi_integral where argument is an array of inputs."""
        methodVals = Fermi_integral(self.args, self.order1)
        testTrue = np.allclose(methodVals, self.Trues, rtol=1e-16, atol=0.0)
        errStr = f"Fermi integral value should be {self.Trues} and not {methodVals}."
        assert testTrue, errStr

    def test_invalid_type(self):
        """
        Test whether `TypeError` is raised when an invalid argument
        type is passed to `~plasmapy.mathematics.Fermi_integral`.
        """
        with pytest.raises(TypeError):
            Fermi_integral([1, 2, 3], self.order1)
