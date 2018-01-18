import pytest
from ..parsing import (
    _get_standard_symbol,
    _is_special_particle,
    _case_insensitive_aliases,
    _case_sensitive_aliases,
    _parse_and_check_atomic_input,)

aliases_and_symbols = [
    ('electron', 'e-'),
    ('beta-', 'e-'),
    ('beta+', 'e+'),
    ('positron', 'e+'),
    ('proton', 'p'),
    ('', ''),
    (5, 5),
    ('deuterium+', 'D 1+'),
    ('deuterium 1+', 'D 1+'),
    ('tritium +1', 'T 1+'),
    ('alpha', 'He-4 2+'),
    ('D+', 'D 1+'),
    ('Deuterium', 'D'),
    ('deuteron', 'D 1+'),
    ('triton', 'T 1+'),
    ('muon', 'mu-'),
    ('antimuon', 'mu+'),
    ('tau particle', 'tau-'),
    ('antitau', 'tau+'),
    ('p+', 'p'),
]


@pytest.mark.parametrize("alias,symbol", aliases_and_symbols)
def test_get_standard_symbol(alias, symbol):
    """Test that _get_standard_symbol correctly takes in aliases and
    returns the corresponding symbols, and returns the original argument
    if the argument does not correspond to an alias."""
    result = _get_standard_symbol(alias)
    assert result == symbol, \
        (f"_get_standard_symbol({alias}) returns {result}, which differs "
         f"from the expected symbol of {symbol}.\n\n"
         f"_case_insensitive_aliases:\n{_case_insensitive_aliases}\n\n"
         f"_case_sensitive_aliases:\n{_case_sensitive_aliases}")


alias_dictionaries = [_case_sensitive_aliases, _case_insensitive_aliases]


@pytest.mark.parametrize("alias_dict", alias_dictionaries)
def test_alias_dict_properties(alias_dict):
    """Test properties of the alias dictionaries."""

    for key in alias_dict.keys():
        assert isinstance(key, str), \
            (f"The following key should be a string, but isn't: {key}\n\n"
             f"The entire dictionary is:\n\n{alias_dict}")

    for value in alias_dict.values():
        assert isinstance(value, str), \
            (f"The following value should be a string, but isn't: {value}\n\n"
             f"The entire dictionary is:\n\n{alias_dict}")



# (arg, kwargs, expected)
parse_check_table = [

    ('He', {'Z': 1, 'mass_numb': 4},
     {'symbol': 'He-4 1+',
      'element': 'He',
      'isotope': 'He-4',
      'ion': 'He-4 1+',
      'mass_numb': 4,
      'Z': 1}),

    (1, {},
     {'symbol': 'H',
      'element': 'H',
      'isotope': None,
      'ion': None,
      'Z': None,
      'mass_numb': None}),

    ('H', {'mass_numb': 2},
     {'symbol': 'D',
      'element': 'H',
      'isotope': 'D',
      'ion': None,
      'Z': None,
      'mass_numb': 2}),

    (2, {},
     {'symbol': 'He',
      'element': 'He',
      'isotope': None,
      'ion': None,
      'Z': None,
      'mass_numb': None}),
]


@pytest.mark.parametrize('arg, kwargs, expected', parse_check_table)
def test_parse_and_check_atomic_input(arg, kwargs, expected):
    result = _parse_and_check_atomic_input(arg, **kwargs)

    assert result == expected, (
        "Error in _parse_and_check_atomic_input.\n"
        "The resulting dictionary is:\n\n"
        f"{result}\n\n"
        "whereas the expected dictionary is:\n\n"
        f"{expected}\n"
    )
