"""
Tests for 'check` decorators (i.e. decorators that only check values but do not
change them).
"""
import inspect
import numpy as np
import pytest

from astropy import units as u
from plasmapy.constants import c
from plasmapy.utils.decorators.checks import (
    _check_quantity,
    _check_relativistic,
    check_values,
    check_quantity,
    check_relativistic,
    CheckValues,
)
from plasmapy.utils.exceptions import (RelativityWarning, RelativityError)
from unittest import mock


# ----------------------------------------------------------------------------------------
# Test Decorator class `CheckValues` and decorator `check_values`
# ----------------------------------------------------------------------------------------
class TestCheckValues:
    """
    Tests for decorator :func:`~plasmapy.utils.decorators.checks.check_values` and
    decorator class :class:`~plasmapy.utils.decorators.checks.CheckValues`.
    """

    @staticmethod
    def foo(x, y):
        return x + y

    def test_cv_default_check_values(self):
        """Test the default check dictionary for CheckValues"""
        cv = CheckValues()
        assert hasattr(cv, '_check_item_defaults')
        assert isinstance(cv._check_item_defaults, dict)
        _defaults = [('can_be_negative', True),
                     ('can_be_complex', False),
                     ('can_be_inf', True),
                     ('can_be_nan', True),
                     ('none_shall_pass', False)]
        for key, val in _defaults:
            assert cv._check_item_defaults[key] == val

    def test_cv_method__check_value(self):
        """
        Test functionality/behavior of the `_check_value` method on `CheckValues`.
        This method does the actual checking of the argument values and should be
        called by `CheckValues.__call__()`.
        """
        # setup wrapped function
        cv = CheckValues()
        wfoo = cv(self.foo)

        assert hasattr(cv, '_check_value')

        # -- Test 'can_be_negative' check --
        check = cv._check_item_defaults
        args = [-5,
                -5.0,
                np.array([-1, 2]),
                np.array([-3., 2.]),
                -3 * u.cm,
                np.array([-4., 3.]) * u.kg]
        for arg in args:
            # can_be_negative == False
            check['can_be_negative'] = False
            with pytest.raises(ValueError):
                cv._check_value(arg, 'arg', **check)

            # can_be_negative == True
            check['can_be_negative'] = True
            assert cv._check_value(arg, 'arg', **check) is None

        # -- Test 'can_be_complex' check --
        check = cv._check_item_defaults
        args = [complex(5),
                complex(2, 3),
                np.complex(3.),
                complex(4., 2.) * u.cm,
                np.array([complex(4, 5), complex(1)]) * u.kg]
        for arg in args:
            # can_be_complex == False
            check['can_be_complex'] = False
            with pytest.raises(ValueError):
                cv._check_value(arg, 'arg', **check)

            # can_be_negative == True
            check['can_be_complex'] = True
            assert cv._check_value(arg, 'arg', **check) is None

        # -- Test 'can_be_inf' check --
        check = cv._check_item_defaults
        args = [np.inf,
                np.inf * u.cm,
                np.array([1., 2., np.inf, 10.]),
                np.array([1., 2., np.inf, np.inf]) * u.kg]
        for arg in args:
            # can_be_inf == False
            check['can_be_inf'] = False
            with pytest.raises(ValueError):
                cv._check_value(arg, 'arg', **check)

            # can_be_inf == True
            check['can_be_inf'] = True
            assert cv._check_value(arg, 'arg', **check) is None

        # -- Test 'can_be_inf' check --
        check = cv._check_item_defaults
        args = [np.nan,
                np.nan * u.cm,
                np.array([1., 2., np.nan, 10.]),
                np.array([1., 2., np.nan, np.nan]) * u.kg]
        for arg in args:
            # can_be_nan == False
            check['can_be_nan'] = False
            with pytest.raises(ValueError):
                cv._check_value(arg, 'arg', **check)

            # can_be_nan == True
            check['can_be_nan'] = True
            assert cv._check_value(arg, 'arg', **check) is None

        # -- Test 'none_shall_pass' check --
        check = cv._check_item_defaults
        args = [None]
        for arg in args:
            # none_shall_pass == False
            check['none_shall_pass'] = False
            with pytest.raises(ValueError):
                cv._check_value(arg, 'arg', **check)

            # none_shall_pass == True
            check['none_shall_pass'] = True
            assert cv._check_value(arg, 'arg', **check) is None

    def test_cv_preserves_signature(self):
        """Test CheckValues preserves signature of wrapped function."""
        # I'd like to directly dest the @preserve_signature is used (??)

        wfoo = CheckValues()(self.foo)
        assert hasattr(wfoo, '__signature__')
        assert wfoo.__signature__ == inspect.signature(self.foo)

    @mock.patch(CheckValues.__module__ + '.' + CheckValues.__qualname__,
                side_effect=CheckValues, autospec=True)
    def test_decorator_definition(self, mock_cv_class):
        """
        Test that :func:`~plasmapy.utils.decorators.checks.check_values` is
        properly defined.
        """
        # create mock function (mock_foo) from function to mock (self.foo)
        mock_foo = mock.Mock(side_effect=self.foo, name='mock_foo', autospec=True)
        mock_foo.__name__ = 'mock_foo'
        mock_foo.__signature__ = inspect.signature(self.foo)

        # define various possible check dicts
        check_list = [
            {'none_shall_pass': False},
            {'can_be_negative': True,
             'can_be_nan': False},
            {'can_be_negative': False,
             'can_be_complex': True,
             'can_be_inf': False,
             'can_be_nan': False,
             'none_shall_pass': True}
        ]

        # Examine various checks and their pass-through
        for check in check_list:
            # -- decorate like a traditional function --
            # decorate
            wfoo = check_values(mock_foo, **{'x': check, 'y': check})

            # tests
            assert wfoo(2, 3) == 5
            assert mock_cv_class.called
            assert mock_foo.called
            assert len(mock_cv_class.call_args) == 2
            assert mock_cv_class.call_args[0] == ()
            assert isinstance(mock_cv_class.call_args[1], dict)
            assert sorted(mock_cv_class.call_args[1].keys()) == ['x', 'y']
            assert mock_cv_class.call_args[1]['x'] == check
            assert mock_cv_class.call_args[1]['y'] == check

            # reset
            mock_cv_class.reset_mock()
            mock_foo.reset_mock()

            # -- decorate like "sugar" syntax --
            #
            #  @check_values(x=check)
            #      def foo(x):
            #          return x
            #
            # decorate
            wfoo = check_values(**{'x': check, 'y': check})(mock_foo)

            # tests
            assert wfoo(2, 3) == 5
            assert mock_cv_class.called
            assert mock_foo.called
            assert len(mock_cv_class.call_args) == 2
            assert mock_cv_class.call_args[0] == ()
            assert isinstance(mock_cv_class.call_args[1], dict)
            assert sorted(mock_cv_class.call_args[1].keys()) == ['x', 'y']
            assert mock_cv_class.call_args[1]['x'] == check
            assert mock_cv_class.call_args[1]['y'] == check

            # reset
            mock_cv_class.reset_mock()
            mock_foo.reset_mock()

        # -- decorate like "sugar" syntax w/o checks--
        #
        #  @check_values
        #      def foo(x):
        #          return x
        #
        # decorate
        wfoo = check_values(mock_foo)

        # tests
        assert wfoo(2, 3) == 5
        assert mock_cv_class.called
        assert mock_foo.called
        assert len(mock_cv_class.call_args) == 2
        assert mock_cv_class.call_args[0] == ()
        assert mock_cv_class.call_args[1] == {}

        # reset
        mock_cv_class.reset_mock()
        mock_foo.reset_mock()


# ----------------------------------------------------------------------------------------
# Test Decorator `check_quantity` (& function `_check_quantity`
# ----------------------------------------------------------------------------------------
# (value, units, error)
quantity_error_examples_default = [
    # exceptions associated with the units keyword
    (5 * u.T, 5 * u.T, TypeError),
    (5 * u.T, 5, TypeError),
    (5 * u.T, [u.T, 1], TypeError),
    (5 * u.T, [1, u.m], TypeError),
    (u.T, u.J, TypeError),
    (3 * u.m / u.s, u.m, u.UnitConversionError),
    (5j * u.K, u.K, ValueError),
]

# (value, units, can_be_negative, can_be_complex, can_be_inf, error)
quantity_error_examples_non_default = [
    (-5 * u.K, u.K, False, False, True, ValueError),
    (np.inf * u.K, u.K, True, False, False, ValueError)
]

# (value, units)
quantity_valid_examples_default = [
    # check basic functionality
    (5 * u.T, u.T),
    (3 * u.m / u.s, u.m / u.s),
    (3 * u.m / u.s, [u.m / u.s]),
    (3 * u.m / u.s ** 2, [u.m / u.s, u.m / (u.s ** 2)]),
    (3 * u.km / u.yr, u.m / u.s),
    # check temperature in units of energy per particle (e.g., eV)
    (5 * u.eV, u.K),
    (5 * u.K, u.eV),
    (5 * u.keV, [u.m, u.K]),
    # check keywords relating to numerical values
    (np.inf * u.T, u.T)
]

# (value, units, can_be_negative, can_be_complex, can_be_inf)
quantity_valid_examples_non_default = [
    (5j * u.m, u.m, True, True, True)
]


# Tests for _check_quantity
@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf, error",
    quantity_error_examples_non_default)
def test__check_quantity_errors_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf, error):
    with pytest.raises(error):
        _check_quantity(value, 'arg', 'funcname', units,
                        can_be_negative=can_be_negative,
                        can_be_complex=can_be_complex,
                        can_be_inf=can_be_inf)


def test__check_quantity_warns_on_casting():
    with pytest.warns(u.UnitsWarning):
        _check_quantity(5, 'arg', 'funcname', u.m,)

@pytest.mark.parametrize(
    "value, units, error", quantity_error_examples_default)
def test__check_quantity_errors_default(value, units, error):
    with pytest.raises(error):
        _check_quantity(value, 'arg', 'funcname', units)


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf",
    quantity_valid_examples_non_default)
def test__check_quantity_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf):
    _check_quantity(value, 'arg', 'funcname', units,
                    can_be_negative=can_be_negative,
                    can_be_complex=can_be_complex,
                    can_be_inf=can_be_inf)


@pytest.mark.parametrize("value, units", quantity_valid_examples_default)
def test__check_quantity_default(value, units):
    _check_quantity(value, 'arg', 'funcname', units)


# Tests for check_quantity decorator
@pytest.mark.parametrize(
    "value, units, error", quantity_error_examples_default)
def test_check_quantity_decorator_errors_default(value, units, error):

    @check_quantity(x={"units": units})
    def func(x):
        return x

    with pytest.raises(error):
        func(value)


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf, error",
    quantity_error_examples_non_default)
def test_check_quantity_decorator_errors_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf, error):

    @check_quantity(
       x={"units": units, "can_be_negative": can_be_negative,
          "can_be_complex": can_be_complex, "can_be_inf": can_be_inf}
    )
    def func(x):
        return x

    with pytest.raises(error):
        func(value)


@pytest.mark.parametrize("value, units", quantity_valid_examples_default)
def test_check_quantity_decorator_default(value, units):

    @check_quantity(x={"units": units})
    def func(x):
        return x

    func(value)


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf",
    quantity_valid_examples_non_default)
def test_check_quantity_decorator_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf):

    @check_quantity(x={"units": units, "can_be_negative": can_be_negative,
                       "can_be_complex": can_be_complex, "can_be_inf": can_be_inf}
                    )
    def func(x):
        return x

    func(value)


def test_check_quantity_decorator_missing_validated_params():

    @check_quantity(
        x={"units": u.m},
        y={"units": u.s}
    )
    def func(x):
        return x

    with pytest.raises(TypeError) as e:
        func(1 * u.m)

    assert "Call to func is missing validated params y" == str(e.value)


def test_check_quantity_decorator_two_args_default():

    @check_quantity(
        x={"units": u.m},
        y={"units": u.s}
    )
    def func(x, y):
        return x / y

    func(1 * u.m, 1 * u.s)


def test_check_quantity_decorator_two_args_not_default():

    @check_quantity(
        x={"units": u.m, "can_be_negative": False},
        y={"units": u.s}
    )
    def func(x, y):
        return x / y

    with pytest.raises(ValueError):
        func(-1 * u.m, 2 * u.s)


def test_check_quantity_decorator_two_args_one_kwargs_default():

    @check_quantity(
        x={"units": u.m},
        y={"units": u.s},
        z={"units": u.eV}
    )
    def func(x, y, another, z=10 * u.eV):
        return x * y * z

    func(1 * u.m, 1 * u.s, 10 * u.T)


def test_check_quantity_decorator_two_args_one_kwargs_not_default():

    @check_quantity(
        x={"units": u.m},
        y={"units": u.s, "can_be_negative": False},
        z={"units": u.eV, "can_be_inf": False}
    )
    def func(x, y, z=10 * u.eV):
        return x * y * z

    with pytest.raises(ValueError):
        func(1 * u.m, 1 * u.s, z=np.inf * u.eV)


class TestCheckQuantityNoneShallPass:
    @check_quantity(x={"units": u.m, "none_shall_pass": True})
    def func(self, x=None):
        if x is None:
            return 0 * u.m
        return x

    def test_incorrect_units(self):
        with pytest.raises(u.UnitConversionError):
            self.func(1*u.s)

    def test_none_to_zero(self):
        assert self.func(None) == 0*u.m


# ----------------------------------------------------------------------------------------
# Test Decorator `check_relativistic` (& function `_check_relativistic`
# ----------------------------------------------------------------------------------------
# (speed, betafrac)
non_relativistic_speed_examples = [
    (0 * u.m / u.s, 0.1),
    (0.0099999 * c, 0.1),
    (-0.009 * c, 0.1),
    (5 * u.AA / u.Gyr, 0.1)
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
    _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac, error", relativistic_error_examples)
def test__check_relativistic_errors(speed, betafrac, error):
    with pytest.raises(error):
        _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac",
                         relativistic_warning_examples)
def test__check_relativistic_warnings(speed, betafrac):
    with pytest.warns(RelativityWarning):
        _check_relativistic(speed, 'f', betafrac=betafrac)


# Tests for check_relativistic decorator
@pytest.mark.parametrize("speed, betafrac", non_relativistic_speed_examples)
def test_check_relativistic_decorator(speed, betafrac):

    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize(
    "speed",
    [item[0] for item in non_relativistic_speed_examples])
def test_check_relativistic_decorator_no_args(speed):

    @check_relativistic
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize(
    "speed",
    [item[0] for item in non_relativistic_speed_examples])
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
