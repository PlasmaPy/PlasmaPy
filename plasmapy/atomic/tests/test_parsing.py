import pytest
from ..parsing import _parse_and_check_atomic_input


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
