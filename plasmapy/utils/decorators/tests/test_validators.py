"""
Tests for 'validate` decorators (i.e. decorators that check objects and change them
when possible).
"""
import inspect
import pytest

from astropy import units as u
from functools import cached_property
from typing import Any, Dict, List
from unittest import mock

from plasmapy.utils.decorators.checks import CheckUnits, CheckValues
from plasmapy.utils.decorators.validators import (
    validate_class_attributes,
    validate_quantities,
    ValidateQuantities,
)


# ----------------------------------------------------------------------------------------
# Test Decorator class `ValidateQuantities` and decorator `validate_quantities`
# ----------------------------------------------------------------------------------------
class TestValidateQuantities:
    """
    Test for decorator
    :class:`~plasmapy.utils.decorators.validators.validate_quantities` and decorator
    class :class:`~plasmapy.utils.decorators.validators.ValidateQuantities`.
    """

    unit_check_defaults = CheckUnits._CheckUnits__check_defaults  # type: Dict[str, Any]
    value_check_defaults = (
        CheckValues._CheckValues__check_defaults
    )  # type: Dict[str, Any]
    check_defaults = {**unit_check_defaults, **value_check_defaults}

    @staticmethod
    def foo(x):
        return x

    @staticmethod
    def foo_anno(x: u.cm):
        return x

    def test_inheritance(self):
        assert issubclass(ValidateQuantities, CheckUnits)
        assert issubclass(ValidateQuantities, CheckValues)

    def test_vq_method__get_validations(self):
        # method must exist
        assert hasattr(ValidateQuantities, "validations")
        assert hasattr(ValidateQuantities, "_get_validations")

        # setup default validations
        default_validations = {
            **self.check_defaults.copy(),
            "units": [self.check_defaults["units"]],
        }

        # setup test cases
        # 'setup' = arguments for `_get_validations`
        # 'output' = expected return from `_get_validations`
        # 'raises' = if `_get_validations` raises an Exception
        # 'warns' = if `_get_validations` issues a warning
        #
        _cases = [
            {
                "descr": "typical call...using 'can_be_negative'",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {"x": {"units": u.cm, "can_be_negative": False}},
                },
                "output": {"x": {"units": [u.cm], "can_be_negative": False}},
            },
            {
                "descr": "typical call...using 'none_shall_pass'",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {"x": {"units": u.cm, "none_shall_pass": True}},
                },
                "output": {"x": {"units": [u.cm], "none_shall_pass": True}},
            },
            {
                "descr": "call w/o value validations",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {"x": {"units": u.cm}},
                },
                "output": {"x": {"units": [u.cm]}},
            },
            {
                "descr": "call w/o unit validations",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {"x": {"can_be_inf": False}},
                },
                "raises": ValueError,
            },
            {
                "descr": "'none_shall_pass' defined w/ validations",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {"x": {"units": [u.cm, None]}},
                },
                "output": {"x": {"units": [u.cm], "none_shall_pass": True}},
            },
            {
                "descr": "units are defined via function annotations",
                "setup": {
                    "function": self.foo_anno,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {},
                },
                "output": {"x": {"units": [u.cm]}},
            },
            {
                "descr": "define 'validations_on_return",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {"validations_on_return": {"units": [u.cm, None]}},
                },
                "output": {
                    "validations_on_return": {"units": [u.cm], "none_shall_pass": True}
                },
            },
            {
                "descr": "'none_shall_pass' is inconsistently doubly defined'",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {
                        "x": {"units": [u.cm, None], "none_shall_pass": False}
                    },
                },
                "raises": ValueError,
            },
            {
                "descr": "define both validations on args and validations_on_return",
                "setup": {
                    "function": self.foo,
                    "args": (5,),
                    "kwargs": {},
                    "validations": {
                        "x": {"units": [u.cm], "none_shall_pass": False},
                        "validations_on_return": {
                            "units": [u.cm],
                            "can_be_zero": False,
                        },
                    },
                },
                "output": {
                    "x": {"units": [u.cm], "none_shall_pass": False},
                    "validations_on_return": {"units": [u.cm], "can_be_zero": False},
                },
            },
        ]  # type: List[Dict[str, Any]]

        for case in _cases:
            sig = inspect.signature(case["setup"]["function"])
            args = case["setup"]["args"]
            kwargs = case["setup"]["kwargs"]
            bound_args = sig.bind(*args, **kwargs)

            vq = ValidateQuantities(**case["setup"]["validations"])
            vq.f = case["setup"]["function"]
            if "warns" in case:
                with pytest.warns(case["warns"]):
                    validations = vq._get_validations(bound_args)
            elif "raises" in case:
                with pytest.raises(case["raises"]):
                    vq._get_validations(bound_args)
                continue
            else:
                validations = vq._get_validations(bound_args)

            # only expected argument validations exist
            assert sorted(validations.keys()) == sorted(case["output"].keys())

            # if validation key-value not specified then default is assumed
            for arg_name in case["output"].keys():
                arg_validations = validations[arg_name]

                for key in default_validations.keys():
                    if key in case["output"][arg_name]:
                        val = case["output"][arg_name][key]
                    else:
                        val = default_validations[key]

                    assert arg_validations[key] == val

        # method calls `_get_unit_checks` and `_get_value_checks`
        with mock.patch.object(
            CheckUnits, "_get_unit_checks", return_value={}
        ) as mock_cu_get, mock.patch.object(
            CheckValues, "_get_value_checks", return_value={}
        ) as mock_cv_get:
            vq = ValidateQuantities(x=u.cm)
            vq.f = self.foo
            sig = inspect.signature(self.foo)
            bound_args = sig.bind(5)

            assert vq._get_validations(bound_args) == {}
            assert mock_cu_get.called
            assert mock_cv_get.called

    def test_vq_method__validate_quantity(self):

        # method must exist
        assert hasattr(ValidateQuantities, "_validate_quantity")

        # setup default validations
        default_validations = {
            **self.check_defaults.copy(),
            "units": [self.check_defaults["units"]],
        }

        # setup test cases
        # 'setup' = arguments for `_get_validations`
        # 'output' = expected return from `_get_validations`
        # 'raises' = if `_get_validations` raises an Exception
        # 'warns' = if `_get_validations` issues a warning
        #
        _cases = [
            # typical call
            {
                "input": {
                    "args": (5 * u.cm, "arg"),
                    "validations": {**default_validations, "units": [u.cm]},
                },
                "output": 5 * u.cm,
            },
            # argument does not have units, but only one is specified
            {
                "input": {
                    "args": (5, "arg"),
                    "validations": {**default_validations, "units": [u.cm]},
                },
                "output": 5 * u.cm,
                "warns": u.UnitsWarning,
            },
            # argument does not have units and multiple unit validations specified
            {
                "input": {
                    "args": (5, "arg"),
                    "validations": {**default_validations, "units": [u.cm, u.km]},
                },
                "raises": TypeError,
            },
            # units can NOT be applied to argument
            {
                "input": {
                    "args": ({}, "arg"),
                    "validations": {**default_validations, "units": [u.cm]},
                },
                "raises": TypeError,
            },
            # argument has a standard unit conversion
            {
                "input": {
                    "args": (5.0 * u.cm, "arg"),
                    "validations": {**default_validations, "units": [u.km]},
                },
                "output": (5.0 * u.cm).to(u.km),
            },
            # return value is None and not allowed
            {
                "input": {
                    "args": (None, "validations_on_return"),
                    "validations": {
                        **default_validations,
                        "units": [u.cm],
                        "none_shall_pass": False,
                    },
                },
                "raises": ValueError,
            },
            # 'pass_equivalent_units' is True and unit conversion is not performed
            {
                "input": {
                    "args": (5 * u.cm, "arg"),
                    "validations": {
                        **default_validations,
                        "units": [u.km],
                        "pass_equivalent_units": True,
                    },
                },
                "output": 5 * u.cm,
            },
        ]

        # setup wrapped function
        vq = ValidateQuantities()
        vq.f = self.foo

        # perform tests
        for ii, case in enumerate(_cases):
            arg, arg_name = case["input"]["args"]
            validations = case["input"]["validations"]

            if "warns" in case:
                with pytest.warns(case["warns"]):
                    _result = vq._validate_quantity(arg, arg_name, validations)
            elif "raises" in case:
                with pytest.raises(case["raises"]):
                    vq._validate_quantity(arg, arg_name, validations)
                continue
            else:
                _result = vq._validate_quantity(arg, arg_name, validations)

            assert _result == case["output"]

        # method calls `_check_unit_core` and `_check_value`
        case = {
            "input": (5.0 * u.cm, u.cm, {**default_validations, "units": [u.cm]}),
            "output": 5.0 * u.cm,
        }
        with mock.patch.object(
            CheckUnits, "_check_unit_core", return_value=(5 * u.cm, u.cm, None, None)
        ) as mock_cu_checks, mock.patch.object(
            CheckValues, "_check_value", return_value=None
        ) as mock_cv_checks:

            args = case["input"][0:2]
            validations = case["input"][2]

            vq = ValidateQuantities(**validations)
            vq.f = self.foo

            assert vq._validate_quantity(*args, validations) == case["output"]
            assert mock_cu_checks.called
            assert mock_cv_checks.called

    def test_vq_preserves_signature(self):
        """Test `ValidateQuantities` preserves signature of wrapped function."""
        # I'd like to directly test the @preserve_signature is used (??)

        wfoo = ValidateQuantities()(self.foo_anno)
        assert hasattr(wfoo, "__signature__")
        assert wfoo.__signature__ == inspect.signature(self.foo_anno)

    def test_vq_called_as_decorator(self):
        """
        Test behavior of `ValidateQuantities.__call__` (i.e. used as a decorator).
        """
        # setup test cases
        # 'setup' = arguments for `CheckUnits` and wrapped function
        # 'output' = expected return from wrapped function
        # 'raises' = if an Exception is expected to be raised
        # 'warns' = if a warning is expected to be issued
        #
        _cases = [
            {
                "descr": "clean execution",
                "setup": {
                    "function": self.foo,
                    "args": (2 * u.cm,),
                    "kwargs": {},
                    "validations": {"x": u.cm, "validations_on_return": u.cm},
                },
                "output": 2 * u.cm,
            },
            {
                "descr": "call with unit conversion",
                "setup": {
                    "function": self.foo,
                    "args": (2 * u.cm,),
                    "kwargs": {},
                    "validations": {"x": u.cm, "validations_on_return": u.um},
                },
                "output": (2 * u.cm).to(u.um),
            },
            {
                "descr": "argument fails checks",
                "setup": {
                    "function": self.foo,
                    "args": (2 * u.cm,),
                    "kwargs": {},
                    "validations": {"x": u.g, "validations_on_return": u.cm},
                },
                "raises": u.UnitTypeError,
            },
            {
                "descr": "return fails checks",
                "setup": {
                    "function": self.foo,
                    "args": (2 * u.cm,),
                    "kwargs": {},
                    "validations": {"x": u.cm, "validations_on_return": u.kg},
                },
                "raises": u.UnitTypeError,
            },
            {
                "descr": "decomposed units are still converted",
                "setup": {
                    "function": self.foo,
                    "args": (2 * u.kg * u.m / u.s**2,),
                    "kwargs": {},
                    "validations": {"x": u.N},
                },
                "output": (2 * u.kg * u.m / u.s**2).to(u.N),
                "extra assert": lambda x: x.unit.to_string() == u.N.to_string(),
            },
        ]

        # perform tests
        for case in _cases:
            validations = case["setup"]["validations"]
            func = case["setup"]["function"]
            wfoo = ValidateQuantities(**validations)(func)

            args = case["setup"]["args"]
            kwargs = case["setup"]["kwargs"]

            if "raises" in case:
                with pytest.raises(case["raises"]):
                    wfoo(*args, **kwargs)
                continue

            _result = wfoo(*args, **kwargs)
            assert _result == case["output"]
            if "extra assert" in case:
                assert case["extra assert"](_result)

        # __call__ calls methods `_get_validations` and `_validate_quantity`
        #
        # setup default validations
        default_validations = self.check_defaults.copy()
        default_validations["units"] = [default_validations.pop("units")]
        default_validations["equivalencies"] = [
            default_validations.pop("equivalencies")
        ]
        validations = {
            "x": {**default_validations, "units": [u.cm]},
            "validations_on_return": {**default_validations, "units": [u.cm]},
        }
        with mock.patch.object(
            ValidateQuantities, "_get_validations", return_value=validations
        ) as mock_vq_get, mock.patch.object(
            ValidateQuantities, "_validate_quantity", return_value=5 * u.cm
        ) as mock_vq_validate:

            wfoo = ValidateQuantities(**validations)(self.foo)
            assert wfoo(5 * u.cm) == 5 * u.cm
            assert mock_vq_get.call_count == 1
            assert mock_vq_validate.call_count == len(validations)

        # validation 'checks_on_return' not allowed
        with pytest.raises(TypeError):
            ValidateQuantities(checks_on_return=u.cm)

        # test on class method
        class Foo:
            @ValidateQuantities()
            def __init__(self, y: u.cm):
                self.y = y

            @ValidateQuantities(validations_on_return={"can_be_negative": False})
            def bar(self, x: u.cm) -> u.m:
                return x + self.y

        foo = Foo(-10 * u.cm)
        assert foo.bar(20.0 * u.cm) == 0.1 * u.m
        with pytest.raises(ValueError):
            foo.bar(5 * u.cm)

    @mock.patch(
        ValidateQuantities.__module__ + "." + ValidateQuantities.__qualname__,
        side_effect=ValidateQuantities,
        autospec=True,
    )
    def test_decorator_func_def(self, mock_vq_class):
        """
        Test that :func:`~plasmapy.utils.decorators.validators.validate_quantities` is
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
                    "args": (5 * u.cm,),
                    "kwargs": {},
                    "validations": {"x": u.cm},
                },
                "output": 5 * u.cm,
            },
            # argument and return checks
            {
                "setup": {
                    "args": (-3.0 * u.cm,),
                    "kwargs": {},
                    "validations": {
                        "x": {"units": u.cm, "can_be_negative": True},
                        "validations_on_return": {
                            "units": u.cm,
                            "can_be_negative": True,
                        },
                    },
                },
                "output": -3 * u.cm,
            },
        ]
        for case in _cases:
            for ii in range(2):
                # decorate
                if ii == 0:
                    # functional decorator call
                    wfoo = validate_quantities(mock_foo, **case["setup"]["validations"])
                elif ii == 1:
                    # sugar decorator call
                    #
                    #  @validate_quantities(x=check)
                    #      def foo(x):
                    #          return x
                    #
                    wfoo = validate_quantities(**case["setup"]["validations"])(mock_foo)
                else:
                    continue

                # test
                args = case["setup"]["args"]
                kwargs = case["setup"]["kwargs"]
                assert wfoo(*args, **kwargs) == case["output"]

                assert mock_vq_class.called
                assert mock_foo.called

                assert mock_vq_class.call_args[0] == ()
                assert sorted(mock_vq_class.call_args[1].keys()) == sorted(
                    case["setup"]["validations"].keys()
                )

                for arg_name, validations in case["setup"]["validations"].items():
                    assert mock_vq_class.call_args[1][arg_name] == validations

                # reset
                mock_vq_class.reset_mock()
                mock_foo.reset_mock()


class TestValidateClassAttributes:
    class SampleCase:
        def __init__(self, x: int = None, y: int = None, z: int = None):
            self.x = x
            self.y = y
            self.z = z

        @cached_property
        @validate_class_attributes(expected_attributes=["x"])
        def require_x(self):
            return 0

        @cached_property
        @validate_class_attributes(expected_attributes=["x", "y"])
        def require_x_and_y(self):
            return 0

        @cached_property
        @validate_class_attributes(both_or_either_attributes=[("x", "y")])
        def require_x_or_y(self):
            return 0

        @cached_property
        @validate_class_attributes(
            expected_attributes=["x"], both_or_either_attributes=[("y", "z")]
        )
        def require_x_and_either_y_or_z(self):
            return 0

        @cached_property
        @validate_class_attributes(mutually_exclusive_attributes=[("x", "y")])
        def require_only_either_x_or_y(self):
            return 0

    @pytest.mark.parametrize(
        "test_case_constructor_keyword_arguments",
        [
            {"x": 0},
            {"y": 0},
            {"z": 0},
            {"x": 0, "y": 0},
            {"x": 0, "z": 0},
            {"y": 0, "z": 0},
        ],
    )
    def test_method_errors(self, test_case_constructor_keyword_arguments):
        """
        Test errors raised by the validate_class_attributes decorator.
        """

        test_case = self.SampleCase(**test_case_constructor_keyword_arguments)

        has_x = "x" in test_case_constructor_keyword_arguments.keys()
        has_y = "y" in test_case_constructor_keyword_arguments.keys()
        has_z = "z" in test_case_constructor_keyword_arguments.keys()

        method_return_dictionary = {
            "require_x": has_x,
            "require_x_and_y": has_x and has_y,
            "require_x_or_y": has_x or has_y,
            "require_x_and_either_y_or_z": has_x and (has_y or has_z),
            "require_only_either_x_or_y": (has_x and not has_y)
            or (has_y and not has_x),
        }

        for method, expected_to_return in method_return_dictionary.items():
            if expected_to_return:
                getattr(test_case, method)
            else:
                with pytest.raises(ValueError):
                    getattr(test_case, method)
