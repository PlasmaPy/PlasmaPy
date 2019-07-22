"""
Tests for 'validate` decorators (i.e. decorators that check objects and change them
when possible).
"""
import inspect
import pytest

from astropy import units as u
from plasmapy.utils.decorators.checks import (CheckUnits, CheckValues)
from plasmapy.utils.decorators.validators import (validate_quantities,
                                                  ValidateQuantities)
from typing import (Any, Dict)
from unittest import mock


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
    value_check_defaults = CheckValues._CheckValues__check_defaults  # type: Dict[str, Any]
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
        assert hasattr(ValidateQuantities, '_get_validations')

        # setup default validations
        default_validations = self.check_defaults.copy()
        default_validations['units'] = [default_validations.pop('units')]
        default_validations['equivalencies'] = [default_validations.pop('equivalencies')]

        # setup test cases
        # 'setup' = arguments for `_get_validations`
        # 'output' = expected return from `_get_validations`
        # 'raises' = if `_get_validations` raises an Exception
        # 'warns' = if `_get_validations` issues a warning
        #
        _cases = [
            # typical call
            {'setup': {'function': self.foo,
                       'args': (5, ),
                       'kwargs': {},
                       'validations': {'x': {'units': u.cm, 'can_be_negative': False}},
                       },
             'output': {'x': {'units': [u.cm],
                              'can_be_negative': False}},
             },
            {'setup': {'function': self.foo,
                       'args': (5,),
                       'kwargs': {},
                       'validations': {'x': {'units': u.cm, 'none_shall_pass': True}},
                       },
             'output': {'x': {'units': [u.cm],
                              'none_shall_pass': True}},
             },

            # call w/o value validations
            {'setup': {'function': self.foo,
                       'args': (5,),
                       'kwargs': {},
                       'validations': {'x': {'units': u.cm}},
                       },
             'output': {'x': {'units': [u.cm]}},
             },

            # call w/o unit validations
            {'setup': {'function': self.foo,
                       'args': (5,),
                       'kwargs': {},
                       'validations': {'x': {'can_be_inf': False}},
                       },
             'raises': ValueError,
             },

            # 'none_shall_pass' defined w/ validations
            {'setup': {'function': self.foo,
                       'args': (5,),
                       'kwargs': {},
                       'validations': {'x': {'units': [u.cm, None]}},
                       },
             'output': {'x': {'units': [u.cm],
                              'none_shall_pass': True}},
             },

            # units are defined via function annotations
            {'setup': {'function': self.foo_anno,
                       'args': (5,),
                       'kwargs': {},
                       'validations': {},
                       },
             'output': {'x': {'units': [u.cm]}},
             },

            # define 'validations_on_return'
            {'setup': {'function': self.foo,
                       'args': (5,),
                       'kwargs': {},
                       'validations': {'validations_on_return': {'units': [u.cm, None]}},
                       },
             'output': {'validations_on_return': {'units': [u.cm],
                                                  'none_shall_pass': True}},
             },
        ]

        for case in _cases:
            sig = inspect.signature(case['setup']['function'])
            args = case['setup']['args']
            kwargs = case['setup']['kwargs']
            bound_args = sig.bind(*args, **kwargs)

            vq = ValidateQuantities(**case['setup']['validations'])
            vq.f = case['setup']['function']
            if 'warns' in case:
                with pytest.warns(case['warns']):
                    validations = vq._get_validations(bound_args)
            elif 'raises' in case:
                with pytest.raises(case['raises']):
                    vq._get_validations(bound_args)
                continue
            else:
                validations = vq._get_validations(bound_args)

            # only expected argument validations exist
            assert sorted(validations.keys()) == sorted(case['output'].keys())

            # if validation key-value not specified then default is assumed
            for arg_name in case['output'].keys():
                arg_validations = validations[arg_name]

                for key in default_validations.keys():
                    if key in case['output'][arg_name]:
                        val = case['output'][arg_name][key]
                    else:
                        val = default_validations[key]

                    assert arg_validations[key] == val

        # method calls `_get_unit_checks` and `_get_value_checks
        with mock.patch.object(CheckUnits, '_get_unit_checks', return_value={}) \
                as mock_cu_get, \
                mock.patch.object(CheckValues, '_get_value_checks', return_value={}) \
                        as mock_cv_get:
            vq = ValidateQuantities(x=u.cm)
            vq.f = self.foo
            sig = inspect.signature(self.foo)
            bound_args = sig.bind(5)

            assert vq._get_validations(bound_args) == {}
            assert mock_cu_get.called
            assert mock_cv_get.called
