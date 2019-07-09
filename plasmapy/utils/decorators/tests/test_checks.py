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
    check_units,
    check_values,
    check_quantity,
    check_relativistic,
    CheckUnits,
    CheckValues,
)
from plasmapy.utils.exceptions import (PlasmaPyWarning,
                                       RelativityWarning,
                                       RelativityError)
from typing import (Any, Dict)
from unittest import mock


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
    def foo(x, y):
        return x + y

    def test_cu_default_check_values(self):
        """Test the default check dictionary for CheckUnits."""
        cu = CheckUnits()
        assert hasattr(cu, '_CheckUnits__check_defaults')
        assert isinstance(cu._CheckUnits__check_defaults, dict)
        _defaults = [('units', None),
                     ('equivalencies', None),
                     ('pass_equivalent_units', False)]
        for key, val in _defaults:
            assert cu._CheckUnits__check_defaults[key] == val

    def test_cu_method__check_unit(self):
        """
        Test functionality/behavior of the `_check_unit` method on `CheckUnits`.
        This method does the actual checking of the argument values and should be
        called by `CheckUnits.__call__()`.
        """
        # setup wrapped function
        cu = CheckUnits()
        wfoo = cu(self.foo)

        assert hasattr(cu, '_check_unit')

        # -- Test 'units' check --
        # setup argument checks
        check = self.check_defaults.copy()
        check['units'] = [u.cm]
        check['equivalencies'] = [None]

        # argument does not have units
        with pytest.raises(TypeError):
            cu._check_unit(5., 'arg', **check)

        # argument does not match desired units
        with pytest.raises(u.UnitTypeError):
            cu._check_unit(5. * u.kg, 'arg', **check)

        # argument has equivalent but not matching units
        with pytest.raises(u.UnitTypeError):
            cu._check_unit(5. * u.km, 'arg', **check)

        # argument is equivalent to many specified units but exactly matches one
        check['units'] = [u.cm, u.km]
        check['equivalencies'] = [None] * 2
        assert cu._check_unit(5. * u.km, 'arg', **check) == (5. * u.km, u.km, None)

        # argument is equivalent to many specified units and does NOT exactly matches one
        with pytest.raises(u.UnitTypeError):
            cu._check_unit(5. * u.m, 'arg', **check)

        # -- Test 'none_shall_pass' check --
        # setup argument checks
        check = self.check_defaults.copy()
        check['units'] = [u.cm]
        check['equivalencies'] = [None]

        # argument is None and none_shall_pass = False
        check['none_shall_pass'] = False
        with pytest.raises(ValueError):
            cu._check_unit(None, 'arg', **check)

        # argument is None and none_shall_pass = True
        check['none_shall_pass'] = True
        assert cu._check_unit(None, 'arg', **check) == (None, None, None)

        # -- Test 'pass_equivalent_units' check --
        # setup argument checks
        check = self.check_defaults.copy()
        check['units'] = [u.cm]
        check['equivalencies'] = [None]
        check['pass_equivalent_units'] = True
        assert cu._check_unit(2. * u.km, 'arg', **check) == (2. * u.km, None, None)

    @mock.patch.object(CheckUnits, '_check_unit')
    def test_cu_called_as_decorator(self, mock_cu):
        """
        Test behavior of `CheckUnits` and argument pass-through to
        `CheckUnits._check_unit` when `CheckUnits` is used as a decorator.
        """
        # What is tested?
        # 1.  basic decorator call where...
        #     * arg x units are specified via decorator keyword argument
        #     * arg y units are specified via function annotations (no additional
        #       checks conditions are specified via decorator keyword argument)
        # 2.  basic decorator call where...
        #     * arg x units are specified via decorator keyword argument
        #     * arg y units are specified via function annotations but additional
        #       check conditions are specified via decorator keyword argument
        # 3.  basic decorator call with None values allowed...
        #     * arg x None values are allows via decorator keyword argument
        #     * arg y None values are allowed via function annotations
        #       (i.e. default value)
        #
        # create mock function (mock_foo) from function to mock (foo)
        def foo(x: u.Quantity, y: u.cm):
            return x.value + y.value
        mock_foo = mock.Mock(side_effect=foo, name='mock_foo', autospec=True)
        mock_foo.__name__ = 'mock_foo'
        mock_foo.__signature__ = inspect.signature(foo)

        # default unit check values
        default_checks = self.check_defaults

        # -- basic wrap and call                                                    -- (1)
        # * units for x will be defined by keyword argument via CheckUnits
        # * units for y will be defined by annotations
        #
        # create wrapped function
        xchecks = {'units': [u.cm],
                   'equivalencies': [u.temperature_energy()]}
        cu = CheckUnits(x=xchecks)
        wfoo = cu(mock_foo)

        # basic tests
        assert wfoo(2 * u.cm, 3 * u.cm) == 5
        assert mock_foo.called
        assert mock_cu.called

        # ensure `_check_value` method was called with correct arguments
        assert mock_cu.call_count == 2
        for ii, (arg, arg_name, checks) in enumerate(zip([2 * u.cm, 3 * u.cm],
                                                         ['x', 'y'],
                                                         [xchecks, {'units': [u.cm]}])):
            # test passed arguments
            assert mock_cu.mock_calls[ii][1] == (arg, arg_name)

            # test passed keywords
            for key, val in default_checks.items():
                if key in checks:
                    # if key defined in checks then value should be passed
                    val = checks[key]
                else:
                    # if key NOT defined in checks then default value should be passed
                    val = default_checks[key]

                if not isinstance(val, (bool, list)):
                    val = [val]

                assert mock_cu.mock_calls[ii][2][key] == val

        # reset mocks
        mock_cu.reset_mock()
        mock_foo.reset_mock()

        # -- Argument checks specified, but relying on function annotations         -- (2)
        # -- to define units                                                        --
        # * units for x will be defined by keyword argument via CheckUnits
        # * units for y will be defined by annotations, but other check
        #   conditions will be defined by keyword arguement
        #
        # create wrapped function
        xchecks = {'units': [u.cm],
                   'equivalencies': [u.temperature_energy()]}
        ychecks = {'pass_equivalent_units': False}
        cu = CheckUnits(x=xchecks, y=ychecks)
        wfoo = cu(mock_foo)

        # basic tests
        assert wfoo(2 * u.cm, 3 * u.cm) == 5
        assert mock_foo.called
        assert mock_cu.called

        # ensure `_check_value` method was called with correct arguments
        ychecks['units'] = [u.cm]
        assert mock_cu.call_count == 2
        for ii, (arg, arg_name, checks) in enumerate(zip([2 * u.cm, 3 * u.cm],
                                                         ['x', 'y'],
                                                         [xchecks, ychecks])):
            # test passed arguments
            assert mock_cu.mock_calls[ii][1] == (arg, arg_name)

            # test passed keywords
            for key, val in default_checks.items():
                if key in checks:
                    # if key defined in checks then value should be passed
                    val = checks[key]
                else:
                    # if key NOT defined in checks then default value should be passed
                    val = default_checks[key]

                if not isinstance(val, (bool, list)):
                    val = [val]

                assert mock_cu.mock_calls[ii][2][key] == val

        # reset mocks
        mock_cu.reset_mock()
        mock_foo.reset_mock()

        # -- Decorator calls with None values                                       -- (3)
        #     * arg x None values are allows via decorator keyword argument
        #     * arg y None values are allowed via function annotations
        #       (i.e. default value)
        #
        # create mock function (mock_foo) from function to mock (foo)
        def foo(x: u.Quantity, y: u.cm = None):
            return x.value + y.value
        mock_foo = mock.Mock(side_effect=foo, name='mock_foo', autospec=True)
        mock_foo.__name__ = 'mock_foo'
        mock_foo.__signature__ = inspect.signature(foo)

        # create wrapped function
        xchecks = {'units': [u.cm, None],
                   'equivalencies': [u.temperature_energy()]}
        cu = CheckUnits(x=xchecks)
        wfoo = cu(mock_foo)

        # basic tests
        assert wfoo(2 * u.cm, 3 * u.cm) == 5
        assert mock_foo.called
        assert mock_cu.called

        # ensure `_check_value` method was called with correct arguments
        assert mock_cu.call_count == 2
        for ii, (arg, arg_name, checks) in enumerate(zip([2 * u.cm, 3 * u.cm],
                                                         ['x', 'y'],
                                                         [xchecks, {'units': [u.cm]}])):
            # test passed arguments
            assert mock_cu.mock_calls[ii][1] == (arg, arg_name)

            # test definition of 'none_shall_pass'
            assert mock_cu.mock_calls[ii][2]['none_shall_pass'] == True

        # reset mocks
        mock_cu.reset_mock()
        mock_foo.reset_mock()


    def test_cu_preserves_signature(self):
        """Test CheckValues preserves signature of wrapped function."""
        # I'd like to directly dest the @preserve_signature is used (??)

        wfoo = CheckUnits()(self.foo)
        assert hasattr(wfoo, '__signature__')
        assert wfoo.__signature__ == inspect.signature(self.foo)

    @mock.patch(CheckUnits.__module__ + '.' + CheckUnits.__qualname__,
                side_effect=CheckUnits, autospec=True)
    def test_decorator_func_def(self, mock_cu_class):
        """
        Test that :func:`~plasmapy.utils.decorators.checks.check_units` is
        properly defined.
        """
        # create mock function (mock_foo) from function to mock (self.foo)
        mock_foo = mock.Mock(side_effect=self.foo, name='mock_foo', autospec=True)
        mock_foo.__name__ = 'mock_foo'
        mock_foo.__signature__ = inspect.signature(self.foo)

        # define various possible check dicts
        check_list = [
            u.cm,
            {'units': u.cm},
            {'units': u.cm,
             'equivalencies': u.temperature()},
        ]

        # Examine various checks and their pass-through
        for check in check_list:
            # -- decorate like a traditional function --
            # decorate
            wfoo = check_units(mock_foo, **{'x': check, 'y': check})

            # tests
            assert wfoo(2 * u.cm, 3 * u.cm) == 5 * u.cm
            assert mock_cu_class.called
            assert mock_foo.called
            assert len(mock_cu_class.call_args) == 2
            assert mock_cu_class.call_args[0] == ()
            assert isinstance(mock_cu_class.call_args[1], dict)
            assert sorted(mock_cu_class.call_args[1].keys()) == ['x', 'y']
            assert mock_cu_class.call_args[1]['x'] == check
            assert mock_cu_class.call_args[1]['y'] == check

            # reset
            mock_cu_class.reset_mock()
            mock_foo.reset_mock()

            # -- decorate like "sugar" syntax --
            #
            #  @check_units(x=check)
            #      def foo(x):
            #          return x
            #
            # decorate
            wfoo = check_units(**{'x': check, 'y': check})(mock_foo)

            # tests
            assert wfoo(2 * u.cm, 3 * u.cm) == 5 * u.cm
            assert mock_cu_class.called
            assert mock_foo.called
            assert len(mock_cu_class.call_args) == 2
            assert mock_cu_class.call_args[0] == ()
            assert isinstance(mock_cu_class.call_args[1], dict)
            assert sorted(mock_cu_class.call_args[1].keys()) == ['x', 'y']
            assert mock_cu_class.call_args[1]['x'] == check
            assert mock_cu_class.call_args[1]['y'] == check

            # reset
            mock_cu_class.reset_mock()
            mock_foo.reset_mock()

        # -- decorate like "sugar" syntax w/o checks --
        #
        #  @check_values
        #      def foo(x):
        #          return x
        #
        # decorate
        wfoo = check_units(mock_foo)

        # tests
        assert wfoo(2, 3) == 5
        assert mock_cu_class.called
        assert mock_foo.called
        assert len(mock_cu_class.call_args) == 2
        assert mock_cu_class.call_args[0] == ()
        assert mock_cu_class.call_args[1] == {}

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

    @staticmethod
    def foo(x, y):
        return x + y

    def test_cv_default_check_values(self):
        """Test the default check dictionary for CheckValues"""
        cv = CheckValues()
        assert hasattr(cv, '_check_defaults')
        assert isinstance(cv._check_defaults, dict)
        _defaults = [('can_be_negative', True),
                     ('can_be_complex', False),
                     ('can_be_inf', True),
                     ('can_be_nan', True),
                     ('none_shall_pass', False)]
        for key, val in _defaults:
            assert cv._check_defaults[key] == val

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
        check = cv._check_defaults.copy()
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
        check = cv._check_defaults.copy()
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
        check = cv._check_defaults.copy()
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
        check = cv._check_defaults.copy()
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
        check = cv._check_defaults.copy()
        args = [None]
        for arg in args:
            # none_shall_pass == False
            check['none_shall_pass'] = False
            with pytest.raises(ValueError):
                cv._check_value(arg, 'arg', **check)

            # none_shall_pass == True
            check['none_shall_pass'] = True
            assert cv._check_value(arg, 'arg', **check) is None

    @mock.patch.object(CheckValues, '_check_value')
    def test_cv_called_as_decorator(self, mock_cv):
        # create mock function (mock_foo) from function to mock (self.foo)
        mock_foo = mock.Mock(side_effect=self.foo, name='mock_foo', autospec=True)
        mock_foo.__name__ = 'mock_foo'
        mock_foo.__signature__ = inspect.signature(self.foo)

        # -- basic wrap and call --
        # create wrapped function
        xchecks = {}
        ychecks = {'none_shall_pass': True,
                   'can_be_negative': False}
        cv = CheckValues(x=xchecks, y=ychecks)
        default_checks = cv._check_defaults.copy()
        wfoo = cv(mock_foo)

        # basic tests
        assert wfoo(2, 3) == 5
        assert mock_foo.called
        assert mock_cv.called

        # ensure `_check_value` method was called with correct arguments
        assert mock_cv.call_count == 2
        for ii, (arg, arg_name, checks) in enumerate(zip([2, 3],
                                                         ['x', 'y'],
                                                         [xchecks, ychecks])):
            # test passed arguments
            assert mock_cv.mock_calls[ii][1] == (arg, arg_name)

            # test passed keywords
            for key, val in default_checks.items():
                if key in checks:
                    # if key defined in checks then value should be passed
                    assert mock_cv.mock_calls[ii][2][key] == checks[key]
                else:
                    # if key NOT defined in checks then default value should be passed
                    assert mock_cv.mock_calls[ii][2][key] == default_checks[key]

        # reset mocks
        mock_cv.reset_mock()
        mock_foo.reset_mock()

        # -- decorator specifies too many arguments --
        # if the number of checked arguments is not less than or equal to the number
        # of function arguments, then a PlasmaPyWarning is raised
        cv = CheckValues(x=xchecks, y=ychecks, z={'none_shall_pass': True})
        with pytest.warns(PlasmaPyWarning):
            wfoo = cv(self.foo)
            assert wfoo(2, 3) == 5

        # reset mocks
        mock_cv.reset_mock()
        mock_foo.reset_mock()

        # -- decorated function has *arg and **kwargs arguments --
        # values passed via *args and **kwargs are ignored by CheckValues
        #
        # create mock function (mock_foo) from function to mock (foo)
        def foo(x, *args, y=3, **kwargs):
            return x + y
        mock_foo = mock.Mock(side_effect=foo, name='mock_foo', autospec=True)
        mock_foo.__name__ = 'mock_foo'
        mock_foo.__signature__ = inspect.signature(foo)

        # create wrapped function
        xchecks = {}
        ychecks = {'none_shall_pass': True,
                   'can_be_negative': False}
        cv = CheckValues(x=xchecks, y=ychecks, z={'none_shall_pass': True})
        wfoo = cv(mock_foo)

        # test
        with pytest.warns(PlasmaPyWarning):
            assert wfoo(5, 'hello', y=3, z=None) == 8
        assert mock_foo.called
        assert mock_cv.call_count == 2

        # reset mocks
        mock_cv.reset_mock()
        mock_foo.reset_mock()

    def test_cv_preserves_signature(self):
        """Test CheckValues preserves signature of wrapped function."""
        # I'd like to directly dest the @preserve_signature is used (??)

        wfoo = CheckValues()(self.foo)
        assert hasattr(wfoo, '__signature__')
        assert wfoo.__signature__ == inspect.signature(self.foo)

    @mock.patch(CheckValues.__module__ + '.' + CheckValues.__qualname__,
                side_effect=CheckValues, autospec=True)
    def test_decorator_func_def(self, mock_cv_class):
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

        # -- decorate like "sugar" syntax w/o checks --
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
