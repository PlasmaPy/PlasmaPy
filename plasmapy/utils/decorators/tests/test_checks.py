"""
Tests for 'check` decorators (i.e. decorators that only check objects but do not
change them).
"""
import inspect
import numpy as np
import pytest

from astropy import units as u
from astropy.constants import c
from plasmapy.utils.decorators.checks import (
    _check_quantity,
    _check_relativistic,
    check_units,
    check_values,
    check_quantity,
    check_relativistic,
    CheckBase,
    CheckUnits,
    CheckValues,
)
from plasmapy.utils.exceptions import (PlasmaPyWarning,
                                       RelativityWarning,
                                       RelativityError)
from typing import (Any, Dict)
from unittest import mock


# ----------------------------------------------------------------------------------------
# Test Decorator class `CheckBase`
# ----------------------------------------------------------------------------------------
class TestCheckBase:
    """
    Test for decorator class :class:`~plasmapy.utils.decorators.checks.CheckBase`.
    """
    def test_for_members(self):
        assert hasattr(CheckUnits, 'checks')

    def test_checks(self):
        _cases = [
            {'input': (None, {'x': 1, 'y': 2}),
             'output': {'x': 1, 'y': 2}},
            {'input': (6, {'x': 1, 'y': 2}),
             'output': {'x': 1, 'y': 2, 'checks_on_return': 6}},
        ]
        for case in _cases:
            cb = CheckBase(checks_on_return=case['input'][0], **case['input'][1])
            assert cb.checks == case['output']


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
    def foo_partial_anno(x: u.Quantity, y: u.cm):
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
        assert hasattr(cu, '_CheckUnits__check_defaults')
        assert isinstance(cu._CheckUnits__check_defaults, dict)
        _defaults = [('units', None),
                     ('equivalencies', None),
                     ('pass_equivalent_units', False),
                     ('none_shall_pass', False)]
        for key, val in _defaults:
            assert cu._CheckUnits__check_defaults[key] == val

    def test_cu_method__get_unit_checks(self):
        """
        Test functionality/behavior of the method `_get_unit_checks` on `CheckUnits`.
        This method reviews the decorator `checks` arguments and wrapped function
        annotations to build a complete checks dictionary.
        """
        # methods must exist
        assert hasattr(CheckUnits, '_get_unit_checks')

        # setup default checks
        default_checks = self.check_defaults.copy()
        default_checks['units'] = [default_checks.pop('units')]
        default_checks['equivalencies'] = [default_checks.pop('equivalencies')]

        # setup test cases
        # 'setup' = arguments for `_get_unit_checks`
        # 'output' = expected return from `_get_unit_checks`
        # 'raises' = if `_get_unit_checks` raises an Exception
        # 'warns' = if `_get_unit_checks` issues a warning
        #
        equivs = [[u.temperature_energy()]]
        _cases = [
            # x units are defined via decorator kwarg of CheckUnits
            # y units are defined via function annotations, additional checks
            #   thru CheckUnits kwarg
            {'setup': {'function': self.foo_partial_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': {'units': [u.cm],
                                        'equivalencies': equivs[0][0]},
                                  },
                       },
             'output': {'x': {'units': [u.cm],
                              'equivalencies': equivs[0]},
                        'y': {'units': [u.cm]},
                        },
             },

            # x units are defined via decorator kwarg of CheckUnits
            # y units are defined via function annotations, additional checks
            #   thru CheckUnits kwarg
            {'setup': {'function': self.foo_partial_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': {'units': [u.cm],
                                        'equivalencies': equivs[0]},
                                  'y': {'pass_equivalent_units': False}
                                  },
                       },
             'output': {'x': {'units': [u.cm],
                              'equivalencies': equivs[0]},
                        'y': {'units': [u.cm],
                              'pass_equivalent_units': False},
                        },
             },

            # number of checked arguments exceed number of function arguments
            {'setup': {'function': self.foo_partial_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': {'units': [u.cm]},
                                  'y': {'units': [u.cm]},
                                  'z': {'units': [u.cm]},
                                  },
                       },
             'warns': PlasmaPyWarning,
             'output': {'x': {'units': [u.cm]},
                        'y': {'units': [u.cm]},
                        },
             },

            # arguments passed via *args and **kwargs are ignored
            {'setup': {'function': self.foo_stars,
                       'args': (2 * u.cm, 'hello'),
                       'kwargs': {'z': None},
                       'checks': {'x': {'units': [u.cm]},
                                  'y': {'units': [u.cm]},
                                  'z': {'units': [u.cm]},
                                  },
                       },
             'warns': PlasmaPyWarning,
             'output': {'x': {'units': [u.cm]},
                        'y': {'units': [u.cm]},
                        },
             },

            # arguments arguments can be None values
            {'setup': {'function': self.foo_with_none,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': {'units': [u.cm, None]}},
                       },
             'output': {'x': {'units': [u.cm],
                              'none_shall_pass': True},
                        'y': {'units': [u.cm],
                              'none_shall_pass': True},
                        },
             },

            # checks and annotations do not specify units
            {'setup': {'function': self.foo_no_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': {'pass_equivalent_units': True}},
                       },
             'raises': ValueError,
             },

            # checks specify too many equivalency lists
            {'setup': {'function': self.foo_no_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': {'units': u.cm,
                                        'equivalencies': [u.temperature(),
                                                          u.temperature_energy()]}},
                       },
             'raises': ValueError,
             },

            # units are directly assigned to the check kwarg
            {'setup': {'function': self.foo_partial_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': u.cm},
                       },
             'output': {'x': {'units': [u.cm]},
                        'y': {'units': [u.cm]},
                        },
             },

            # return units are assigned via checks
            {'setup': {'function': self.foo_return_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'checks_on_return': u.km},
                       },
             'output': {'checks_on_return': {'units': [u.km]},
                        },
             },

            # return units are assigned via annotations
            {'setup': {'function': self.foo_return_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {},
                       },
             'output': {'checks_on_return': {'units': [u.um]},
                        },
             },

            # return units are not specified but other checks are
            {'setup': {'function': self.foo_no_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'checks_on_return': {'pass_equivalent_units': True}},
                       },
             'raises': ValueError,
             },

        ]

        # perform tests
        for case in _cases:
            sig = inspect.signature(case['setup']['function'])
            bound_args = sig.bind(*case['setup']['args'], **case['setup']['kwargs'])

            cu = CheckUnits(**case['setup']['checks'])
            cu.f = case['setup']['function']
            if 'warns' in case:
                with pytest.warns(case['warns']):
                    checks = cu._get_unit_checks(bound_args)
            elif 'raises' in case:
                with pytest.raises(case['raises']):
                    cu._get_unit_checks(bound_args)
                continue
            else:
                checks = cu._get_unit_checks(bound_args)

            # only expected argument checks exist
            assert sorted(checks.keys()) == sorted(case['output'].keys())

            # if check key-value not specified then default is assumed
            for arg_name in case['output'].keys():
                arg_checks = checks[arg_name]

                for key in default_checks.keys():
                    if key in case['output'][arg_name]:
                        val = case['output'][arg_name][key]
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
        assert hasattr(CheckUnits, '_check_unit')
        assert hasattr(CheckUnits, '_check_unit_core')

        # setup default checks
        check = self.check_defaults.copy()
        check['units'] = [u.cm]
        check['equivalencies'] = [None]

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
            {'input': (5., 'arg', check),
             'output': (None, None, None, TypeError)},

            # argument does match desired units
            # * set arg_name = 'checks_on_return' to cover if-else statement
            #   in initializing error string
            {'input': (5. * u.kg, 'checks_on_return', check),
             'output': (None, None, None, u.UnitTypeError)},

            # argument has equivalent but not matching unit
            {'input': (5. * u.km, 'arg', check),
             'output': (5. * u.km, u.cm, None, u.UnitTypeError)},

            # argument is equivalent to many specified units but exactly matches one
            {'input': (5. * u.km,
                       'arg', {**check,
                               'units': [u.cm, u.km],
                               'equivalencies': [None] * 2}),
             'output': (5. * u.km, u.km, None, None)},

            # argument is equivalent to many specified units and
            # does NOT exactly match one
            {'input': (5. * u.m,
                       'arg', {**check,
                               'units': [u.cm, u.km],
                               'equivalencies': [None] * 2}),
             'output': (None, None, None, u.UnitTypeError)},

            # argument has attr unit but unit does not have is_equivalent
            {'input': (MyQuantity, 'arg', check),
             'output': (None, None, None, TypeError)},
        ]

        # add cases for 'none_shall_pass' checks
        _cases.extend([
            # argument is None and none_shall_pass = False
            {'input': (None, 'arg', {**check, 'none_shall_pass': False}),
             'output': (None, None, None, ValueError)},

            # argument is None and none_shall_pass = True
            {'input': (None, 'arg', {**check, 'none_shall_pass': True}),
             'output': (None, None, None, None)},

        ])

        # add cases for 'pass_equivalent_units' checks
        _cases.extend([
            # argument is equivalent to 1 to unit,
            # does NOT exactly match the unit,
            # and 'pass_equivalent_units' = True and argument
            {'input': (5. * u.km, 'arg', {**check, 'pass_equivalent_units': True}),
             'output': (5. * u.km, u.cm, None, None)},

            # argument is equivalent to more than 1 unit,
            # does NOT exactly match any unit,
            # and 'pass_equivalent_units' = True and argument
            {'input': (5. * u.km, 'arg', {**check,
                                          'units': [u.cm, u.m],
                                          'equivalencies': [None] * 2,
                                          'pass_equivalent_units': True}),
             'output': (5. * u.km, None, None, None)},
        ])

        # setup wrapped function
        cu = CheckUnits()
        cu.f = self.foo_no_anno

        # perform tests
        for case in _cases:
            arg, arg_name, arg_checks = case['input']
            _results = cu._check_unit_core(arg, arg_name, **arg_checks)
            assert _results[0:3] == case['output'][0:3]

            if _results[3] is None:
                assert _results[3] is case['output'][3]
                assert cu._check_unit(arg, arg_name, **arg_checks) is None
            else:
                assert isinstance(_results[3], case['output'][3])
                with pytest.raises(case['output'][3]):
                    cu._check_unit(arg, arg_name, **arg_checks)

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
            {'setup': {'function': self.foo_no_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': u.cm,
                                  'y': u.cm,
                                  'checks_on_return': u.cm}},
             'output': 5 * u.cm,
             },

            # argument fails checks
            {'setup': {'function': self.foo_no_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': u.g,
                                  'y': u.cm,
                                  'checks_on_return': u.cm}},
             'raises': u.UnitTypeError,
             },

            # return fails checks
            {'setup': {'function': self.foo_no_anno,
                       'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': u.cm,
                                  'y': u.cm,
                                  'checks_on_return': u.km}},
             'raises': u.UnitTypeError,
             },
        ]

        # test
        for case in _cases:
            wfoo = CheckUnits(**case['setup']['checks'])(case['setup']['function'])

            args = case['setup']['args']
            kwargs = case['setup']['kwargs']

            if 'raises' in case:
                with pytest.raises(case['raises']):
                    wfoo(*args, **kwargs)
            else:
                assert wfoo(*args, **kwargs) == case['output']

    def test_cu_preserves_signature(self):
        """Test `CheckValues` preserves signature of wrapped function."""
        # I'd like to directly dest the @preserve_signature is used (??)

        wfoo = CheckUnits()(self.foo_no_anno)
        assert hasattr(wfoo, '__signature__')
        assert wfoo.__signature__ == inspect.signature(self.foo_no_anno)

    @mock.patch(CheckUnits.__module__ + '.' + CheckUnits.__qualname__,
                side_effect=CheckUnits, autospec=True)
    def test_decorator_func_def(self, mock_cu_class):
        """
        Test that :func:`~plasmapy.utils.decorators.checks.check_units` is
        properly defined.
        """
        # create mock function (mock_foo) from function to mock (self.foo_no_anno)
        mock_foo = mock.Mock(side_effect=self.foo_no_anno, name='mock_foo',
                             autospec=True)
        mock_foo.__name__ = 'mock_foo'
        mock_foo.__signature__ = inspect.signature(self.foo_no_anno)

        # setup test cases
        # 'setup' = arguments for `check_units` and wrapped function
        # 'output' = expected return from wrapped function
        # 'raises' = a raised Exception is expected
        # 'warns' = an issued warning is expected
        #
        _cases = [
            # only argument checks
            {'setup': {'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': u.cm,
                                  'y': u.cm}
                       },
             'output': 5 * u.cm,
             },

            # argument and return checks
            {'setup': {'args': (2 * u.cm, 3 * u.cm),
                       'kwargs': {},
                       'checks': {'x': u.cm,
                                  'checks_on_return': u.cm}
                       },
             'output': 5 * u.cm,
             },
        ]
        for case in _cases:
            for ii in range(2):
                # decorate
                if ii == 0:
                    # functional decorator call
                    wfoo = check_units(mock_foo, **case['setup']['checks'])
                elif ii == 1:
                    # sugar decorator call
                    #
                    #  @check_units(x=check)
                    #      def foo(x):
                    #          return x
                    #
                    wfoo = check_units(**case['setup']['checks'])(mock_foo)
                else:
                    continue

                # test
                args = case['setup']['args']
                kwargs = case['setup']['kwargs']
                assert wfoo(*args, **kwargs) == case['output']

                assert mock_cu_class.called
                assert mock_foo.called

                assert len(mock_cu_class.call_args) == len(case['setup']['checks'])
                assert mock_cu_class.call_args[0] == ()
                assert (sorted(mock_cu_class.call_args[1].keys())
                        == sorted(case['setup']['checks'].keys()))

                for arg_name, checks in case['setup']['checks'].items():
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
        assert hasattr(cv, '_CheckValues__check_defaults')
        assert isinstance(cv._CheckValues__check_defaults, dict)
        _defaults = [('can_be_negative', True),
                     ('can_be_complex', False),
                     ('can_be_inf', True),
                     ('can_be_nan', True),
                     ('none_shall_pass', False)]
        for key, val in _defaults:
            assert cv._CheckValues__check_defaults[key] == val

    def test_cv_method__get_value_checks(self):
        """
        Test functionality/behavior of the method `_get_value_checks` on `CheckValues`.
        This method reviews the decorator `checks` arguments to build a complete
        checks dictionary.
        """
        # methods must exist
        assert hasattr(CheckValues, '_get_value_checks')

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
            {'setup': {'function': self.foo,
                       'args': (2, 3),
                       'kwargs': {},
                       'checks': {'x': {'can_be_negative': False,
                                        'can_be_complex': True,
                                        'can_be_inf': False},
                                  'checks_on_return': {'can_be_nan': False,
                                                       'none_shall_pass': True},
                                  },
                       },
             'output': {'x': {'can_be_negative': False,
                              'can_be_complex': True,
                              'can_be_inf': False},
                        'checks_on_return': {'can_be_nan': False,
                                             'none_shall_pass': True},
                        },
             },

            # arguments passed via *args and **kwargs are ignored
            {'setup': {'function': self.foo_stars,
                       'args': (2, 'hello'),
                       'kwargs': {'z': None},
                       'checks': {'x': {'can_be_negative': False},
                                  'y': {'can_be_inf': False},
                                  'z': {'none_shall_pass': True},
                                  },
                       },
             'output': {'x': {'can_be_negative': False},
                        'y': {'can_be_inf': False},
                        },
             'warns': PlasmaPyWarning,
             },

            # check argument is not a dictionary (default is assumed)
            {'setup': {'function': self.foo,
                       'args': (2, 3),
                       'kwargs': {},
                       'checks': {'x': u.cm},
                       },
             'output': {'x': {}},
             },
        ]

        # perform tests
        for case in _cases:
            sig = inspect.signature(case['setup']['function'])
            args = case['setup']['args']
            kwargs = case['setup']['kwargs']
            bound_args = sig.bind(*args, **kwargs)

            cv = CheckValues(**case['setup']['checks'])
            cv.f = case['setup']['function']
            if 'warns' in case:
                with pytest.warns(case['warns']):
                    checks = cv._get_value_checks(bound_args)
            elif 'raises' in case:
                with pytest.raises(case['raises']):
                    cv._get_value_checks(bound_args)
                continue
            else:
                checks = cv._get_value_checks(bound_args)

            # only expected keys exist
            assert sorted(checks.keys()) == sorted(case['output'].keys())

            # if check key-value not specified then default is assumed
            for arg_name in case['output'].keys():
                arg_checks = checks[arg_name]

                for key in default_checks.keys():
                    if key in case['output'][arg_name]:
                        val = case['output'][arg_name][key]
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
        assert hasattr(cv, '_check_value')

        # setup default checks
        default_checks = self.check_defaults.copy()

        # setup test cases
        # 'setup' = arguments for `CheckUnits` and wrapped function
        # 'raises' = if an Exception is expected to be raised
        # 'warns' = if a warning is expected to be issued
        #
        _cases = [
            # tests for check 'can_be_negative'
            {'input': {'args': [-5,
                                -5.0,
                                np.array([-1, 2]),
                                np.array([-3., 2.]),
                                -3 * u.cm,
                                np.array([-4., 3.]) * u.kg],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'can_be_negative': False}},
             'raises': ValueError,
             },
            {'input': {'args': [-5,
                                -5.0,
                                np.array([-1, 2]),
                                np.array([-3., 2.]),
                                -3 * u.cm,
                                np.array([-4., 3.]) * u.kg],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'can_be_negative': True}},
             },

            # tests for check 'can_be_complex'
            {'input': {'args': [complex(5),
                                complex(2, 3),
                                np.complex(3.),
                                complex(4., 2.) * u.cm,
                                np.array([complex(4, 5), complex(1)]) * u.kg],
                       'arg_name': 'checks_on_return',
                       'checks': {**default_checks, 'can_be_complex': False}},
             'raises': ValueError,
             },
            {'input': {'args': [complex(5),
                                complex(2, 3),
                                np.complex(3.),
                                complex(4., 2.) * u.cm,
                                np.array([complex(4, 5), complex(1)]) * u.kg],
                       'arg_name': 'checks_on_return',
                       'checks': {**default_checks, 'can_be_complex': True}},
             },

            # tests for check 'can_be_inf'
            {'input': {'args': [np.inf,
                                np.inf * u.cm,
                                np.array([1., 2., np.inf, 10.]),
                                np.array([1., 2., np.inf, np.inf]) * u.kg],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'can_be_inf': False}},
             'raises': ValueError,
             },
            {'input': {'args': [np.inf,
                                np.inf * u.cm,
                                np.array([1., 2., np.inf, 10.]),
                                np.array([1., 2., np.inf, np.inf]) * u.kg],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'can_be_inf': True}},
             },

            # tests for check 'can_be_nan'
            {'input': {'args': [np.nan,
                                np.nan * u.cm,
                                np.array([1., 2., np.nan, 10.]),
                                np.array([1., 2., np.nan, np.nan]) * u.kg],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'can_be_nan': False}},
             'raises': ValueError,
             },
            {'input': {'args': [np.nan,
                                np.nan * u.cm,
                                np.array([1., 2., np.nan, 10.]),
                                np.array([1., 2., np.nan, np.nan]) * u.kg],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'can_be_nan': True}},
             },

            # tests for check 'none_shall_pass'
            {'input': {'args': [None],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'none_shall_pass': False}},
             'raises': ValueError,
             },
            {'input': {'args': [None],
                       'arg_name': 'arg',
                       'checks': {**default_checks, 'none_shall_pass': True}},
             },
        ]

        # test
        for case in _cases:
            arg_name = case['input']['arg_name']
            checks = case['input']['checks']

            for arg in case['input']['args']:
                if 'raises' in case:
                    with pytest.raises(case['raises']):
                        cv._check_value(arg, arg_name, **checks)
                elif 'warns' in case:
                    with pytest.warns(case['warns']):
                        cv._check_value(arg, arg_name, **checks)
                else:
                    assert cv._check_value(arg, arg_name, **checks) is None

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
            {'setup': {'function': self.foo,
                       'args': (2, -3),
                       'kwargs': {},
                       'checks': {'x': {'can_be_negative': True},
                                  'y': {'can_be_negative': True},
                                  'checks_on_return': {'can_be_negative': True}}},
             'output': -1,
             },

            # argument fails checks
            {'setup': {'function': self.foo,
                       'args': (2, -3),
                       'kwargs': {},
                       'checks': {'x': {'can_be_negative': True},
                                  'y': {'can_be_negative': False},
                                  'checks_on_return': {'can_be_negative': True}}},
             'raises': ValueError,
             },

            # return fails checks
            {'setup': {'function': self.foo,
                       'args': (2, -3),
                       'kwargs': {},
                       'checks': {'x': {'can_be_negative': True},
                                  'y': {'can_be_negative': True},
                                  'checks_on_return': {'can_be_negative': False}}},
             'raises': ValueError,
             },
        ]

        # test
        for case in _cases:
            wfoo = CheckValues(**case['setup']['checks'])(case['setup']['function'])

            args = case['setup']['args']
            kwargs = case['setup']['kwargs']

            if 'raises' in case:
                with pytest.raises(case['raises']):
                    wfoo(*args, **kwargs)
            else:
                assert wfoo(*args, **kwargs) == case['output']

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

        # setup test cases
        # 'setup' = arguments for `check_units` and wrapped function
        # 'output' = expected return from wrapped function
        # 'raises' = a raised Exception is expected
        # 'warns' = an issued warning is expected
        #
        _cases = [
            # only argument checks
            {'setup': {'args': (-4, 3),
                       'kwargs': {},
                       'checks': {'x': {'can_be_negative': True},
                                  'y': {'can_be_nan': False}}
                       },
             'output': -1,
             },

            # argument and return checks
            {'setup': {'args': (-4, 3),
                       'kwargs': {},
                       'checks': {'x': {'can_be_negative': True},
                                  'checks_on_return': {'can_be_negative': True}}
                       },
             'output': -1,
             },
        ]
        for case in _cases:
            for ii in range(2):
                # decorate
                if ii == 0:
                    # functional decorator call
                    wfoo = check_values(mock_foo, **case['setup']['checks'])
                elif ii == 1:
                    # sugar decorator call
                    #
                    #  @check_values(x=check)
                    #      def foo(x):
                    #          return x
                    #
                    wfoo = check_values(**case['setup']['checks'])(mock_foo)
                else:
                    continue

                # test
                args = case['setup']['args']
                kwargs = case['setup']['kwargs']
                assert wfoo(*args, **kwargs) == case['output']

                assert mock_cv_class.called
                assert mock_foo.called

                assert len(mock_cv_class.call_args) == len(case['setup']['checks'])
                assert mock_cv_class.call_args[0] == ()
                assert (sorted(mock_cv_class.call_args[1].keys())
                        == sorted(case['setup']['checks'].keys()))

                for arg_name, checks in case['setup']['checks'].items():
                    assert mock_cv_class.call_args[1][arg_name] == checks

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
