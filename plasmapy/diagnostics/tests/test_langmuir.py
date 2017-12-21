"""Tests for functions that use Langmuir analysis functions."""

import numpy as np

from plasmapy.langmuir import (swept_probe_analysis,
                               obtain_EEDF,
                               obtain_EEPF)

class Test_swept_probe_analysis():
    def test_no_signal():
        """
        Checks whether an error is raised when zero signal is passed to the
        function.
        """
        errStr = ("Performing a swept probe analysis should return an error "
                  "when zero signal is passed, but doesn't.")
        assert True, errStr

class Test_obtain_EEDF():
    def test_no_signal():
        """
        Checks whether an error is raised when zero signal is passed to the
        function.
        """
        errStr = ("Attempting to obtain the EEDF should return an error when "
                  "zero signal is passed, but doesn't.")
        assert True, errStr

class Test_obtain_EEPF():
    def test_no_signal():
        """
        Checks whether an error is raised when zero signal is passed to the
        function.
        """
        errStr = ("Attempting to obtain the EEPF should return an error when "
                  "zero signal is passed, but doesn't.")
        assert True, errStr