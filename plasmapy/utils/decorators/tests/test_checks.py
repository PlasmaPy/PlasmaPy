"""
Tests for 'check` decorators (i.e. decorators that only check objects but do not
change them).
"""
import inspect
import numpy as np
import pytest

from astropy import units as u
from astropy.constants import c
from types import LambdaType
from typing import Any, Dict
from unittest import mock

from plasmapy.utils.decorators.checks import (
    _check_relativistic,
    check_relativistic,
    check_units,
    check_values,
    CheckBase,
    CheckUnits,
    CheckValues,
)
from plasmapy.utils.exceptions import (
    PlasmaPyWarning,
    RelativityError,
    RelativityWarning,
)


# ----------------------------------------------------------------------------------------
# Test Decorator class `CheckBase`
# ----------------------------------------------------------------------------------------
class TestCheckBase:
    """
    Test for decorator class :class:`~plasmapy.utils.decorators.checks.CheckBase`.
    """

    def test_for_members(self):
        assert hasattr(CheckUnits, "checks")

    def test_checks(self):
        _cases = [
            {"input": (None, {"x": 1, "y": 2}), "output": {"x": 1, "y": 2}},
            {
                "input": (6, {"x": 1, "y": 2}),
                "output": {"x": 1, "y": 2, "checks_on_return": 6},
            },
        ]
        for case in _cases:
            cb = CheckBase(checks_on_return=case["input"][0], **case["input"][1])
            assert cb.checks == case["output"]


# ----------------------------------------------------------------------------------------
# Test Decorator class `CheckValues` and decorator `check_values`
# ----------------------------------------------------------------------------------------
class TestCheckUnits:
    """
    Tests for decorator :func:`~plasmapy.utils.decorators.checks.check_units` and
    decorator class :class:`~plasmapy.utils.decorators.checks.CheckUnits`.
    """

    check_defaults = CheckUnits._CheckUnits__check_defaults  # type: Dict[str, Any]

    @staticmethod
    def foo_no_anno(x, y):
        return x + y

    @staticmethod
    def foo_partial_anno(x: u.Quantity, y: u.cm) -> u.Quantity:
        return x.value + y.value

    @staticmethod
    def foo_return_anno(x, y) -> u.um:
        return x.value + y.value

    @staticmethod
    def foo_stars(x: u.Quantity, *args, y=3 * u.cm, **kwargs):
        return x.value + y.value

    @staticmethod
    def foo_with_none(x: u.Quantity, y: u.cm = None):
        return x.value + y.value

    def test_inheritance(self):
        assert issubclass(CheckUnits, CheckBase)

    def test_cu_default_check_values(self):
        """Test the default check dictionary for CheckUnits."""
        cu = CheckUnits()
        assert hasattr(cu, "_CheckUnits__check_defaults")
        assert isinstance(cu._CheckUnits__check_defaults, dict)
        _defaults = [
            ("units", None),
            ("equivalencies", None),
            ("pass_equivalent_units", False),
            ("none_shall_pass", False),
        ]
        for key, val in _defaults:
            assert cu._CheckUnits__check_defaults[key] == val

    def test_cu_method__flatten_equivalencies_list(self):
        assert hasattr(CheckUnits, "_flatten_equivalencies_list")

        cu = CheckUnits()
        pairs = [([1, 2, 4], [1, 2, 4]), ([1, 2, (3, 4), [5, 6]], [1, 2, (3, 4), 5, 6])]
        for pair in pairs:
            assert cu._flatten_equivalencies_list(pair[0]) == pair[1]

    def test_cu_method__condition_target_units(self):
        """Test method `CheckUnits._condition_target_units`."""
        assert hasattr(CheckUnits, "_condition_target_units")

        cu = CheckUnits()

        targets = ["cm", u.km, u.Quantity, float]
        conditioned_targets = [u.cm, u.km]
        with pytest.raises(TypeError):
            cu._condition_target_units(targets)

        assert (
            cu._condition_target_units(targets, from_annotations=True)
            == conditioned_targets
        )

        with pytest.raises(ValueError):
            cu._condition_target_units(["five"])

    def test_cu_method__normalize_equivalencies(self):
        """Test method `CheckUnits._normalize_equivalencies`."""
        assert hasattr(CheckUnits, "_normalize_equivalencies")

        cu = CheckUnits()

        assert cu._normalize_equivalencies(None) == []

        # 2 element equivalency
        norme = cu._normalize_equivalencies([(u.cm, u.cm)])
        assert len(norme) == 1
        assert isinstance(norme[0], tuple)
        assert len(norme[0]) == 4
        assert norme[0][0] == norme[0][1]
        assert norme[0][2] == norme[0][3]
        assert isinstance(norme[0][2], LambdaType)
        assert norme[0][1] == norme[0][2](norme[0][0])
        assert norme[0][0] == norme[0][3](norme[0][1])

        # 3 element equivalency
        norme = cu._normalize_equivalencies([(u.cm, u.cm, lambda x: x)])
        assert len(norme) == 1
        assert isinstance(norme[0], tuple)
        assert len(norme[0]) == 4
        assert norme[0][0] == norme[0][1]
        assert norme[0][2] == norme[0][3]
        assert isinstance(norme[0][2], LambdaType)
        assert norme[0][1] == norme[0][2](norme[0][0])
        assert norme[0][0] == norme[0][3](norme[0][1])

        # 3 element equivalency
        norme = cu._normalize_equivalencies(
            [(u.K, u.deg_C, lambda x: x - 273.15, lambda x: x + 273.15)]
        )
        assert len(norme) == 1
        assert isinstance(norme[0], tuple)
        assert len(norme[0]) == 4
        assert norme[0][0] == u.K
        assert norme[0][1] == u.deg_C
        assert isinstance(norme[0][2], LambdaType)
        assert isinstance(norme[0][3], LambdaType)
        for val in [-20.0, 50.0, 195.0]:
            assert norme[0][2](val) == (lambda x: x - 273.15)(val)
            assert norme[0][3](val) == (lambda x: x + 273.15)(val)

        # not a 2, 3, or 4-tuple
        with pytest.raises(ValueError):
            cu._normalize_equivalencies([(u.cm,)])

        # input is not a astropy.unit.Unit
        with pytest.raises(ValueError):
            cu._normalize_equivalencies([("cm", u.cm)])

    def test_cu_method__get_unit_checks(self):
        """
        Test functionality/behavior of the method `_get_unit_checks` on `CheckUnits`.
        This method reviews the decorator `checks` arguments and wrapped function
        annotations to build a complete checks dictionary.
        """
        # methods must exist
        assert hasattr(CheckUnits, "_get_unit_checks")

        # setup default checks
        default_checks = {
            **self.check_defaults.copy(),
            "units": [self.check_defaults["units"]],
        }

        # setup test cases
        # 'setup' = arguments for `_get_unit_checks`
        # 'output' = expected return from `_get_unit_checks`
        # 'raises' = if `_get_unit_checks` raises an Exception
        # 'warns' = if `_get_unit_checks` issues a warning
        #
        equivs = [
            # list of astropy Equivalency objects
            [u.temperature_energy(), u.temperature()],
            # list of equivalencies (pre astropy v3.2.1 style)
            list(u.temperature()),
        ]
        _cases = [
            {
                "descr": "x units are defined via decorator kwarg of CheckUnits\n"
                "y units are defined via decorator annotations, additional\n"
                "  checks thru CheckUnits kwarg",
                "setup": {
                    "function": self.foo_partial_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": {"units": [u.cm], "equivalencies": equivs[0][0]}},
                },
                "output": {
                    "x": {"units": [u.cm], "equivalencies": equivs[0][0]},
                    "y": {"units": [u.cm]},
                },
            },
            {
                "descr": "x units are defined via decorator kwarg of CheckUnits\n"
                "y units are defined via function annotations, additional\n"
                "  checks thru CheckUnits kwarg",
                "setup": {
                    "function": self.foo_partial_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {
                        "x": {"units": [u.cm], "equivalencies": equivs[0]},
                        "y": {"pass_equivalent_units": False},
                    },
                },
                "output": {
                    "x": {
                        "units": [u.cm],
                        "equivalencies": equivs[0][0] + equivs[0][1],
                    },
                    "y": {"units": [u.cm], "pass_equivalent_units": False},
                },
            },
            {
                "descr": "equivalencies are a list instead of astropy Equivalency objects",
                "setup": {
                    "function": self.foo_no_anno,
                    "args": (2 * u.K, 3 * u.K),
                    "kwargs": {},
                    "checks": {
                        "x": {"units": [u.K], "equivalencies": equivs[1][0]},
                        "y": {"units": [u.K], "equivalencies": equivs[1]},
                    },
                },
                "output": {
                    "x": {"units": [u.K], "equivalencies": [equivs[1][0]]},
                    "y": {"units": [u.K], "equivalencies": equivs[1]},
                },
            },
            {
                "descr": "number of checked arguments exceed number of function arguments",
                "setup": {
                    "function": self.foo_partial_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {
                        "x": {"units": [u.cm]},
                        "y": {"units": [u.cm]},
                        "z": {"units": [u.cm]},
                    },
                },
                "warns": PlasmaPyWarning,
                "output": {"x": {"units": [u.cm]}, "y": {"units": [u.cm]}},
            },
            {
                "descr": "arguments passed via *args and **kwargs are ignored",
                "setup": {
                    "function": self.foo_stars,
                    "args": (2 * u.cm, "hello"),
                    "kwargs": {"z": None},
                    "checks": {
                        "x": {"units": [u.cm]},
                        "y": {"units": [u.cm]},
                        "z": {"units": [u.cm]},
                    },
                },
                "warns": PlasmaPyWarning,
                "output": {"x": {"units": [u.cm]}, "y": {"units": [u.cm]}},
            },
            {
                "descr": "arguments can be None values",
                "setup": {
                    "function": self.foo_with_none,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": {"units": [u.cm, None]}},
                },
                "output": {
                    "x": {"units": [u.cm], "none_shall_pass": True},
                    "y": {"units": [u.cm], "none_shall_pass": True},
                },
            },
            {
                "descr": "checks and annotations do not specify units",
                "setup": {
                    "function": self.foo_no_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": {"pass_equivalent_units": True}},
                },
                "raises": ValueError,
            },
            {
                "descr": "units are directly assigned to the check kwarg",
                "setup": {
                    "function": self.foo_partial_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": u.cm},
                },
                "output": {"x": {"units": [u.cm]}, "y": {"units": [u.cm]}},
            },
            {
                "descr": "return units are assigned via checks",
                "setup": {
                    "function": self.foo_no_anno,
                    "args": (2 * u.km, 3 * u.km),
                    "kwargs": {},
                    "checks": {"checks_on_return": u.km},
                },
                "output": {"checks_on_return": {"units": [u.km]}},
            },
            {
                "descr": "return units are assigned via annotations",
                "setup": {
                    "function": self.foo_return_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {},
                },
                "output": {"checks_on_return": {"units": [u.um]}},
            },
            {
                "descr": "return units are assigned via annotations and checks arg, but"
                "are not consistent",
                "setup": {
                    "function": self.foo_return_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"checks_on_return": {"units": u.km}},
                },
                "raises": ValueError,
            },
            {
                "descr": "return units are not specified but other checks are",
                "setup": {
                    "function": self.foo_no_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"checks_on_return": {"pass_equivalent_units": True}},
                },
                "raises": ValueError,
            },
            {
                "descr": "no parameter checks for x are defined, but a non-unit annotation"
                "is used",
                "setup": {
                    "function": self.foo_partial_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {},
                },
                "output": {"y": {"units": [u.cm]}},
            },
            {
                "descr": "parameter checks defined for x but unit checks calculated from"
                "function annotations. Function annotations do NOT define "
                "a proper unit type.",
                "setup": {
                    "function": self.foo_partial_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": {"pass_equivalent_units": True}},
                },
                "raises": ValueError,
            },
            {
                "descr": "parameter checks defined for return argument but unit checks "
                "calculated from function annotations. Function annotations do "
                "NOT define a proper unit type.",
                "setup": {
                    "function": self.foo_partial_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"checks_on_return": {"pass_equivalent_units": True}},
                },
                "raises": ValueError,
            },
        ]

        # perform tests
        for ii, case in enumerate(_cases):
            sig = inspect.signature(case["setup"]["function"])
            bound_args = sig.bind(*case["setup"]["args"], **case["setup"]["kwargs"])

            cu = CheckUnits(**case["setup"]["checks"])
            cu.f = case["setup"]["function"]
            if "warns" in case:
                with pytest.warns(case["warns"]):
                    checks = cu._get_unit_checks(bound_args)
            elif "raises" in case:
                with pytest.raises(case["raises"]):
                    cu._get_unit_checks(bound_args)
                continue
            else:
                checks = cu._get_unit_checks(bound_args)

            # only expected argument checks exist
            assert sorted(checks.keys()) == sorted(case["output"].keys())

            # if check key-value not specified then default is assumed
            for arg_name in case["output"].keys():
                arg_checks = checks[arg_name]

                for key in default_checks.keys():
                    if key in case["output"][arg_name]:
                        val = case["output"][arg_name][key]
                    else:
                        val = default_checks[key]

                    assert arg_checks[key] == val

    def test_cu_method__check_unit(self):
        """
        Test functionality/behavior of the methods `_check_unit` and `_check_unit_core`
        on `CheckUnits`.  These methods do the actual checking of the argument units
        and should be called by `CheckUnits.__call__()`.
        """
        # methods must exist
        assert hasattr(CheckUnits, "_check_unit")
        assert hasattr(CheckUnits, "_check_unit_core")

        # setup default checks
        check = {**self.check_defaults, "units": [u.cm]}
        # check = self.check_defaults.copy()
        # check['units'] = [u.cm]
        # check['equivalencies'] = [None]

        # make a class w/ improper units
        class MyQuantity:
            unit = None

        # setup test cases
        # 'input' = arguments for `_check_unit_core` and `_check_unit`
        # 'output' = expected return from `_check_unit_core`
        #
        # add cases for 'units' checks
        _cases = [
            # argument does not have units
            {"input": (5.0, "arg", check), "output": (None, None, None, TypeError)},
            # argument does match desired units
            # * set arg_name = 'checks_on_return' to cover if-else statement
            #   in initializing error string
            {
                "input": (5.0 * u.kg, "checks_on_return", check),
                "output": (None, None, None, u.UnitTypeError),
            },
            # argument has equivalent but not matching unit
            {
                "input": (5.0 * u.km, "arg", check),
                "output": (5.0 * u.km, u.cm, None, u.UnitTypeError),
            },
            # argument is equivalent to many specified units but exactly matches one
            {
                "input": (5.0 * u.km, "arg", {**check, "units": [u.cm, u.km]}),
                "output": (5.0 * u.km, u.km, None, None),
            },
            # argument is equivalent to many specified units and
            # does NOT exactly match one
            {
                "input": (5.0 * u.m, "arg", {**check, "units": [u.cm, u.km]}),
                "output": (None, None, None, u.UnitTypeError),
            },
            # argument has attr unit but unit does not have is_equivalent
            {
                "input": (MyQuantity, "arg", check),
                "output": (None, None, None, TypeError),
            },
        ]

        # add cases for 'none_shall_pass' checks
        _cases.extend(
            [
                # argument is None and none_shall_pass = False
                {
                    "input": (None, "arg", {**check, "none_shall_pass": False}),
                    "output": (None, None, None, ValueError),
                },
                # argument is None and none_shall_pass = True
                {
                    "input": (None, "arg", {**check, "none_shall_pass": True}),
                    "output": (None, None, None, None),
                },
            ]
        )

        # add cases for 'pass_equivalent_units' checks
        _cases.extend(
            [
                # argument is equivalent to 1 to unit,
                # does NOT exactly match the unit,
                # and 'pass_equivalent_units' = True and argument
                {
                    "input": (
                        5.0 * u.km,
                        "arg",
                        {**check, "pass_equivalent_units": True},
                    ),
                    "output": (5.0 * u.km, u.cm, None, None),
                },
                # argument is equivalent to more than 1 unit,
                # does NOT exactly match any unit,
                # and 'pass_equivalent_units' = True and argument
                {
                    "input": (
                        5.0 * u.km,
                        "arg",
                        {**check, "units": [u.cm, u.m], "pass_equivalent_units": True},
                    ),
                    "output": (5.0 * u.km, None, None, None),
                },
            ]
        )

        # setup wrapped function
        cu = CheckUnits()
        cu.f = self.foo_no_anno

        # perform tests
        for ii, case in enumerate(_cases):
            arg, arg_name, arg_checks = case["input"]
            _results = cu._check_unit_core(arg, arg_name, arg_checks)
            assert _results[0:3] == case["output"][0:3]

            if _results[3] is None:
                assert _results[3] is case["output"][3]
                assert cu._check_unit(arg, arg_name, arg_checks) is None
            else:
                assert isinstance(_results[3], case["output"][3])
                with pytest.raises(case["output"][3]):
                    cu._check_unit(arg, arg_name, arg_checks)

    def test_cu_called_as_decorator(self):
        """
        Test behavior of `CheckUnits.__call__` (i.e. used as a decorator).
        """
        # setup test cases
        # 'setup' = arguments for `CheckUnits` and wrapped function
        # 'output' = expected return from wrapped function
        # 'raises' = if an Exception is expected to be raised
        # 'warns' = if a warning is expected to be issued
        #
        _cases = [
            # clean execution
            {
                "setup": {
                    "function": self.foo_no_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": u.cm, "y": u.cm, "checks_on_return": u.cm},
                },
                "output": 5 * u.cm,
            },
            # argument fails checks
            {
                "setup": {
                    "function": self.foo_no_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": u.g, "y": u.cm, "checks_on_return": u.cm},
                },
                "raises": u.UnitTypeError,
            },
            # return fails checks
            {
                "setup": {
                    "function": self.foo_no_anno,
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": u.cm, "y": u.cm, "checks_on_return": u.km},
                },
                "raises": u.UnitTypeError,
            },
        ]

        # test
        for case in _cases:
            wfoo = CheckUnits(**case["setup"]["checks"])(case["setup"]["function"])

            args = case["setup"]["args"]
            kwargs = case["setup"]["kwargs"]

            if "raises" in case:
                with pytest.raises(case["raises"]):
                    wfoo(*args, **kwargs)
            else:
                assert wfoo(*args, **kwargs) == case["output"]

        # test on class method
        class Foo:
            @CheckUnits()
            def __init__(self, y: u.cm):
                self.y = y

            @CheckUnits(x=u.cm)
            def bar(self, x) -> u.cm:
                return x + self.y

        foo = Foo(10.0 * u.cm)
        assert foo.bar(-3 * u.cm) == 7 * u.cm

    def test_cu_preserves_signature(self):
        """Test `CheckValues` preserves signature of wrapped function."""
        # I'd like to directly test the @preserve_signature is used (??)

        wfoo = CheckUnits()(self.foo_no_anno)
        assert hasattr(wfoo, "__signature__")
        assert wfoo.__signature__ == inspect.signature(self.foo_no_anno)

    @mock.patch(
        CheckUnits.__module__ + "." + CheckUnits.__qualname__,
        side_effect=CheckUnits,
        autospec=True,
    )
    def test_decorator_func_def(self, mock_cu_class):
        """
        Test that :func:`~plasmapy.utils.decorators.checks.check_units` is
        properly defined.
        """
        # create mock function (mock_foo) from function to mock (self.foo_no_anno)
        mock_foo = mock.Mock(
            side_effect=self.foo_no_anno, name="mock_foo", autospec=True
        )
        mock_foo.__name__ = "mock_foo"
        mock_foo.__signature__ = inspect.signature(self.foo_no_anno)

        # setup test cases
        # 'setup' = arguments for `check_units` and wrapped function
        # 'output' = expected return from wrapped function
        # 'raises' = a raised Exception is expected
        # 'warns' = an issued warning is expected
        #
        _cases = [
            # only argument checks
            {
                "setup": {
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": u.cm, "y": u.cm},
                },
                "output": 5 * u.cm,
            },
            # argument and return checks
            {
                "setup": {
                    "args": (2 * u.cm, 3 * u.cm),
                    "kwargs": {},
                    "checks": {"x": u.cm, "checks_on_return": u.cm},
                },
                "output": 5 * u.cm,
            },
        ]
        for case in _cases:
            for ii in range(2):
                # decorate
                if ii == 0:
                    # functional decorator call
                    wfoo = check_units(mock_foo, **case["setup"]["checks"])
                elif ii == 1:
                    # sugar decorator call
                    #
                    #  @check_units(x=check)
                    #      def foo(x):
                    #          return x
                    #
                    wfoo = check_units(**case["setup"]["checks"])(mock_foo)
                else:
                    continue

                # test
                args = case["setup"]["args"]
                kwargs = case["setup"]["kwargs"]
                assert wfoo(*args, **kwargs) == case["output"]

                assert mock_cu_class.called
                assert mock_foo.called

                assert mock_cu_class.call_args[0] == ()
                assert sorted(mock_cu_class.call_args[1].keys()) == sorted(
                    case["setup"]["checks"].keys()
                )

                for arg_name, checks in case["setup"]["checks"].items():
                    assert mock_cu_class.call_args[1][arg_name] == checks

                # reset
                mock_cu_class.reset_mock()
                mock_foo.reset_mock()


# ----------------------------------------------------------------------------------------
# Test Decorator class `CheckValues` and decorator `check_values`
# ----------------------------------------------------------------------------------------
class TestCheckValues:
    """
    Tests for decorator :func:`~plasmapy.utils.decorators.checks.check_values` and
    decorator class :class:`~plasmapy.utils.decorators.checks.CheckValues`.
    """

    check_defaults = CheckValues._CheckValues__check_defaults  # type: Dict[str, bool]

    @staticmethod
    def foo(x, y):
        return x + y

    @staticmethod
    def foo_stars(x, *args, y=3, **kwargs):
        return x + y

    def test_inheritance(self):
        assert issubclass(CheckValues, CheckBase)

    def test_cv_default_check_values(self):
        """Test the default check dictionary for CheckValues"""
        cv = CheckValues()
        assert hasattr(cv, "_CheckValues__check_defaults")
        assert isinstance(cv._CheckValues__check_defaults, dict)
        _defaults = [
            ("can_be_negative", True),
            ("can_be_complex", False),
            ("can_be_inf", True),
            ("can_be_nan", True),
            ("none_shall_pass", False),
            ("can_be_zero", True),
        ]
        for key, val in _defaults:
            assert cv._CheckValues__check_defaults[key] == val

    def test_cv_method__get_value_checks(self):
        """
        Test functionality/behavior of the method `_get_value_checks` on `CheckValues`.
        This method reviews the decorator `checks` arguments to build a complete
        checks dictionary.
        """
        # methods must exist
        assert hasattr(CheckValues, "_get_value_checks")

        # setup default checks
        default_checks = self.check_defaults.copy()

        # setup test cases
        # 'setup' = arguments for `_get_value_checks`
        # 'output' = expected return from `_get_value_checks`
        # 'raises' = if `_get_value_checks` raises an Exception
        # 'warns' = if `_get_value_checks` issues a warning
        #
        _cases = [
            # define some checks
            {
                "setup": {
                    "function": self.foo,
                    "args": (2, 3),
                    "kwargs": {},
                    "checks": {
                        "x": {
                            "can_be_negative": False,
                            "can_be_complex": True,
                            "can_be_inf": False,
                        },
                        "checks_on_return": {
                            "can_be_nan": False,
                            "none_shall_pass": True,
                        },
                    },
                },
                "output": {
                    "x": {
                        "can_be_negative": False,
                        "can_be_complex": True,
                        "can_be_inf": False,
                    },
                    "checks_on_return": {"can_be_nan": False, "none_shall_pass": True},
                },
            },
            # arguments passed via *args and **kwargs are ignored
            {
                "setup": {
                    "function": self.foo_stars,
                    "args": (2, "hello"),
                    "kwargs": {"z": None},
                    "checks": {
                        "x": {"can_be_negative": False},
                        "y": {"can_be_inf": False},
                        "z": {"none_shall_pass": True},
                    },
                },
                "output": {"x": {"can_be_negative": False}, "y": {"can_be_inf": False}},
                "warns": PlasmaPyWarning,
            },
            # check argument is not a dictionary (default is assumed)
            {
                "setup": {
                    "function": self.foo,
                    "args": (2, 3),
                    "kwargs": {},
                    "checks": {"x": u.cm},
                },
                "output": {"x": {}},
            },
        ]

        # perform tests
        for case in _cases:
            sig = inspect.signature(case["setup"]["function"])
            args = case["setup"]["args"]
            kwargs = case["setup"]["kwargs"]
            bound_args = sig.bind(*args, **kwargs)

            cv = CheckValues(**case["setup"]["checks"])
            cv.f = case["setup"]["function"]
            if "warns" in case:
                with pytest.warns(case["warns"]):
                    checks = cv._get_value_checks(bound_args)
            elif "raises" in case:
                with pytest.raises(case["raises"]):
                    cv._get_value_checks(bound_args)
                continue
            else:
                checks = cv._get_value_checks(bound_args)

            # only expected keys exist
            assert sorted(checks.keys()) == sorted(case["output"].keys())

            # if check key-value not specified then default is assumed
            for arg_name in case["output"].keys():
                arg_checks = checks[arg_name]

                for key in default_checks.keys():
                    if key in case["output"][arg_name]:
                        val = case["output"][arg_name][key]
                    else:
                        val = default_checks[key]

                    assert arg_checks[key] == val

    def test_cv_method__check_value(self):
        """
        Test functionality/behavior of the `_check_value` method on `CheckValues`.
        This method does the actual checking of the argument values and should be
        called by `CheckValues.__call__()`.
        """
        # setup wrapped function
        cv = CheckValues()
        wfoo = cv(self.foo)

        # methods must exist
        assert hasattr(cv, "_check_value")

        # setup default checks
        default_checks = self.check_defaults.copy()

        # setup test cases
        # 'setup' = arguments for `CheckUnits` and wrapped function
        # 'raises' = if an Exception is expected to be raised
        # 'warns' = if a warning is expected to be issued
        #
        _cases = [
            # tests for check 'can_be_negative'
            {
                "input": {
                    "args": [
                        -5,
                        -5.0,
                        np.array([-1, 2]),
                        np.array([-3.0, 2.0]),
                        -3 * u.cm,
                        np.array([-4.0, 3.0]) * u.kg,
                    ],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_negative": False},
                },
                "raises": ValueError,
            },
            {
                "input": {
                    "args": [
                        -5,
                        -5.0,
                        np.array([-1, 2]),
                        np.array([-3.0, 2.0]),
                        -3 * u.cm,
                        np.array([-4.0, 3.0]) * u.kg,
                    ],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_negative": True},
                }
            },
            # tests for check 'can_be_complex'
            {
                "input": {
                    "args": [
                        complex(5),
                        complex(2, 3),
                        np.complex128(3.0),
                        complex(4.0, 2.0) * u.cm,
                        np.array([complex(4, 5), complex(1)]) * u.kg,
                    ],
                    "arg_name": "checks_on_return",
                    "checks": {**default_checks, "can_be_complex": False},
                },
                "raises": ValueError,
            },
            {
                "input": {
                    "args": [
                        complex(5),
                        complex(2, 3),
                        np.complex128(3.0),
                        complex(4.0, 2.0) * u.cm,
                        np.array([complex(4, 5), complex(1)]) * u.kg,
                    ],
                    "arg_name": "checks_on_return",
                    "checks": {**default_checks, "can_be_complex": True},
                }
            },
            # tests for check 'can_be_inf'
            {
                "input": {
                    "args": [
                        np.inf,
                        np.inf * u.cm,
                        np.array([1.0, 2.0, np.inf, 10.0]),
                        np.array([1.0, 2.0, np.inf, np.inf]) * u.kg,
                    ],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_inf": False},
                },
                "raises": ValueError,
            },
            {
                "input": {
                    "args": [
                        np.inf,
                        np.inf * u.cm,
                        np.array([1.0, 2.0, np.inf, 10.0]),
                        np.array([1.0, 2.0, np.inf, np.inf]) * u.kg,
                    ],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_inf": True},
                }
            },
            # tests for check 'can_be_nan'
            {
                "input": {
                    "args": [
                        np.nan,
                        np.nan * u.cm,
                        np.array([1.0, 2.0, np.nan, 10.0]),
                        np.array([1.0, 2.0, np.nan, np.nan]) * u.kg,
                    ],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_nan": False},
                },
                "raises": ValueError,
            },
            {
                "input": {
                    "args": [
                        np.nan,
                        np.nan * u.cm,
                        np.array([1.0, 2.0, np.nan, 10.0]),
                        np.array([1.0, 2.0, np.nan, np.nan]) * u.kg,
                    ],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_nan": True},
                }
            },
            # tests for check 'none_shall_pass'
            {
                "input": {
                    "args": [None],
                    "arg_name": "arg",
                    "checks": {**default_checks, "none_shall_pass": False},
                },
                "raises": ValueError,
            },
            {
                "input": {
                    "args": [None],
                    "arg_name": "arg",
                    "checks": {**default_checks, "none_shall_pass": True},
                }
            },
            # tests for check 'can_be_zero'
            {
                "input": {
                    "args": [0, 0 * u.cm, np.arange(0, 3), np.arange(0, 3) * u.kg],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_zero": False},
                },
                "raises": ValueError,
            },
            {
                "input": {
                    "args": [0, 0 * u.cm, np.arange(0, 3), np.arange(0, 3) * u.kg],
                    "arg_name": "arg",
                    "checks": {**default_checks, "can_be_zero": True},
                }
            },
        ]

        # test
        for case in _cases:
            arg_name = case["input"]["arg_name"]
            checks = case["input"]["checks"]

            for arg in case["input"]["args"]:
                if "raises" in case:
                    with pytest.raises(case["raises"]):
                        cv._check_value(arg, arg_name, checks)
                elif "warns" in case:
                    with pytest.warns(case["warns"]):
                        cv._check_value(arg, arg_name, checks)
                else:
                    assert cv._check_value(arg, arg_name, checks) is None

    def test_cv_called_as_decorator(self):
        """
        Test behavior of `CheckValues.__call__` (i.e. used as a decorator).
        """
        # setup test cases
        # 'setup' = arguments for `CheckUnits` and wrapped function
        # 'output' = expected return from wrapped function
        # 'raises' = if an Exception is expected to be raised
        # 'warns' = if a warning is expected to be issued
        #
        _cases = [
            # clean execution
            {
                "setup": {
                    "function": self.foo,
                    "args": (2, -3),
                    "kwargs": {},
                    "checks": {
                        "x": {"can_be_negative": True},
                        "y": {"can_be_negative": True},
                        "checks_on_return": {"can_be_negative": True},
                    },
                },
                "output": -1,
            },
            # argument fails checks
            {
                "setup": {
                    "function": self.foo,
                    "args": (2, -3),
                    "kwargs": {},
                    "checks": {
                        "x": {"can_be_negative": True},
                        "y": {"can_be_negative": False},
                        "checks_on_return": {"can_be_negative": True},
                    },
                },
                "raises": ValueError,
            },
            # return fails checks
            {
                "setup": {
                    "function": self.foo,
                    "args": (2, -3),
                    "kwargs": {},
                    "checks": {
                        "x": {"can_be_negative": True},
                        "y": {"can_be_negative": True},
                        "checks_on_return": {"can_be_negative": False},
                    },
                },
                "raises": ValueError,
            },
        ]

        # test on function
        for case in _cases:
            wfoo = CheckValues(**case["setup"]["checks"])(case["setup"]["function"])

            args = case["setup"]["args"]
            kwargs = case["setup"]["kwargs"]

            if "raises" in case:
                with pytest.raises(case["raises"]):
                    wfoo(*args, **kwargs)
            else:
                assert wfoo(*args, **kwargs) == case["output"]

        # test on class method
        class Foo:
            @CheckValues(y={"can_be_negative": True})
            def __init__(self, y):
                self.y = y

            @CheckValues(
                x={"can_be_negative": True}, checks_on_return={"can_be_negative": False}
            )
            def bar(self, x):
                return x + self.y

        foo = Foo(-5)
        assert foo.bar(6) == 1
        with pytest.raises(ValueError):
            foo.bar(1)

    def test_cv_preserves_signature(self):
        """Test CheckValues preserves signature of wrapped function."""
        # I'd like to directly test the @preserve_signature is used (??)

        wfoo = CheckValues()(self.foo)
        assert hasattr(wfoo, "__signature__")
        assert wfoo.__signature__ == inspect.signature(self.foo)

    @mock.patch(
        CheckValues.__module__ + "." + CheckValues.__qualname__,
        side_effect=CheckValues,
        autospec=True,
    )
    def test_decorator_func_def(self, mock_cv_class):
        """
        Test that :func:`~plasmapy.utils.decorators.checks.check_values` is
        properly defined.
        """
        # create mock function (mock_foo) from function to mock (self.foo)
        mock_foo = mock.Mock(side_effect=self.foo, name="mock_foo", autospec=True)
        mock_foo.__name__ = "mock_foo"
        mock_foo.__signature__ = inspect.signature(self.foo)

        # setup test cases
        # 'setup' = arguments for `check_units` and wrapped function
        # 'output' = expected return from wrapped function
        # 'raises' = a raised Exception is expected
        # 'warns' = an issued warning is expected
        #
        _cases = [
            # only argument checks
            {
                "setup": {
                    "args": (-4, 3),
                    "kwargs": {},
                    "checks": {
                        "x": {"can_be_negative": True},
                        "y": {"can_be_nan": False},
                    },
                },
                "output": -1,
            },
            # argument and return checks
            {
                "setup": {
                    "args": (-4, 3),
                    "kwargs": {},
                    "checks": {
                        "x": {"can_be_negative": True},
                        "checks_on_return": {"can_be_negative": True},
                    },
                },
                "output": -1,
            },
        ]
        for case in _cases:
            for ii in range(2):
                # decorate
                if ii == 0:
                    # functional decorator call
                    wfoo = check_values(mock_foo, **case["setup"]["checks"])
                elif ii == 1:
                    # sugar decorator call
                    #
                    #  @check_values(x=check)
                    #      def foo(x):
                    #          return x
                    #
                    wfoo = check_values(**case["setup"]["checks"])(mock_foo)
                else:
                    continue

                # test
                args = case["setup"]["args"]
                kwargs = case["setup"]["kwargs"]
                assert wfoo(*args, **kwargs) == case["output"]

                assert mock_cv_class.called
                assert mock_foo.called

                assert mock_cv_class.call_args[0] == ()
                assert sorted(mock_cv_class.call_args[1].keys()) == sorted(
                    case["setup"]["checks"].keys()
                )

                for arg_name, checks in case["setup"]["checks"].items():
                    assert mock_cv_class.call_args[1][arg_name] == checks

                # reset
                mock_cv_class.reset_mock()
                mock_foo.reset_mock()


# ----------------------------------------------------------------------------------------
# Test Decorator `check_relativistic` (& function `_check_relativistic`
# ----------------------------------------------------------------------------------------
# (speed, betafrac)
non_relativistic_speed_examples = [
    (0 * u.m / u.s, 0.1),
    (0.0099999 * c, 0.1),
    (-0.009 * c, 0.1),
    (5 * u.AA / u.Gyr, 0.1),
]

# (speed, betafrac, error)
relativistic_error_examples = [
    (u.m / u.s, 0.1, TypeError),
    (51513.35, 0.1, TypeError),
    (5 * u.m, 0.1, u.UnitConversionError),
    (1.0 * c, 0.1, RelativityError),
    (1.1 * c, 0.1, RelativityError),
    (np.inf * u.cm / u.s, 0.1, RelativityError),
    (-1.0 * c, 0.1, RelativityError),
    (-1.1 * c, 0.1, RelativityError),
    (-np.inf * u.cm / u.s, 0.1, RelativityError),
]

# (speed, betafrac, warning)
relativistic_warning_examples = [
    (0.11 * c, 0.1),
    (-0.11 * c, 0.1),
    (2997924581 * u.cm / u.s, 0.1),
    (0.02 * c, 0.01),
]


# Tests for _check_relativistic
@pytest.mark.parametrize("speed, betafrac", non_relativistic_speed_examples)
def test__check_relativisitc_valid(speed, betafrac):
    _check_relativistic(speed, "f", betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac, error", relativistic_error_examples)
def test__check_relativistic_errors(speed, betafrac, error):
    with pytest.raises(error):
        _check_relativistic(speed, "f", betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac", relativistic_warning_examples)
def test__check_relativistic_warnings(speed, betafrac):
    with pytest.warns(RelativityWarning):
        _check_relativistic(speed, "f", betafrac=betafrac)


# Tests for check_relativistic decorator
@pytest.mark.parametrize("speed, betafrac", non_relativistic_speed_examples)
def test_check_relativistic_decorator(speed, betafrac):
    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize("speed", [item[0] for item in non_relativistic_speed_examples])
def test_check_relativistic_decorator_no_args(speed):
    @check_relativistic
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize("speed", [item[0] for item in non_relativistic_speed_examples])
def test_check_relativistic_decorator_no_args_parentheses(speed):
    @check_relativistic()
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize("speed, betafrac, error", relativistic_error_examples)
def test_check_relativistic_decorator_errors(speed, betafrac, error):
    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    with pytest.raises(error):
        speed_func()
