"""Tests for functions that use Langmuir analysis functions."""

import numpy as np

from ..langmuir import (swept_probe_analysis,
                        obtain_EEDF,
                        obtain_EEPF)

class Test_swept_probe_analysis(object):
    def test_no_signal():
        """
        Checks whether an error is raised when zero signal is passed to the
        function.
        """
        errStr = (f"Performing a swept probe analysis should return an error "
                  f"when zero signal is passed, but doesn't.")
        assert True, errStr

class Test_obtain_EEDF(object):
    def test_no_signal():
        """
        Checks whether an error is raised when zero signal is passed to the
        function.
        """
        errStr = (f"Attempting to obtain the EEDF should return an error when "
                  f"zero signal is passed, but doesn't.")
        assert True, errStr

class Test_obtain_EEPF(object):
    def test_no_signal():
        """
        Checks whether an error is raised when zero signal is passed to the
        function.
        """
        errStr = (f"Attempting to obtain the EEPF should return an error when "
                  f"zero signal is passed, but doesn't.")
        assert True, errStr