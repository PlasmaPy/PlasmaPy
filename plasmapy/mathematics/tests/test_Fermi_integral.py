"""Tests the Fermi_integral function in the mathematics module."""

# This file contains experimental usage of unicode characters.

import numpy as np
import pytest

from .. mathematics import Fermi_integral

class Test_Fermi_integral(object):
    def setup_method(self):
        """initializing parameters for tests """
        self.arg1 = 3.889780
        self.True1 = (6.272518847136373-8.673617379884035e-19j)
        self.argFail1 = 3.889781
        self.order1 = 0.5
    def test_known1(self):
        """
        Tests Fermi_integral for expected value.
        """
        methodVal = Fermi_integral(self.arg1, self.order1)
        testTrue = np.isclose(methodVal,
                              self.True1,
                              rtol=1e-16,
                              atol=0.0)
        errStr = (f"Fermi integral value should be {self.True1} and not "
                  f"{methodVal}.")
        assert testTrue, errStr
    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 + 1e-15
        methodVal = Fermi_integral(self.arg1, self.order1)
        testTrue = not np.isclose(methodVal,
                                  fail1,
                                  rtol=1e-16,
                                  atol=0.0)
        errStr = (f"Fermi integral value test gives {methodVal} and should "
                  f"be equal to {fail1}.")
        assert testTrue, errStr
    def test_polog_fail(self):
        """
        Tests whether Fermi_integral() fails due to polylog from mpmath
        not having an implementation for larger argument values.
        """
        with pytest.raises(NotImplementedError):
            Fermi_integral(self.argFail1, self.order1)